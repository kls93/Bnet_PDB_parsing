
# Code run in terminal to update PDB_file_properties, and to identify the
# nearest resolution structures (minimum 1000) about each entry in
# PDB_file_properties_filtered + calculate relevant statistics

# Resolution bin calculation
import os
import shutil
import numpy as np
import pandas as pd
from collections import OrderedDict
from scipy import stats

if os.path.isdir('Resolution_filter'):
    shutil.rmtree('Resolution_filter')
os.mkdir('Resolution_filter')

df = pd.read_csv('PDB_file_properties_filtered.csv')
df = df.reset_index(drop=True)

df_len = df.shape[0]
sub_df_stats = OrderedDict({'PDB code': [np.nan]*df_len,
                            'Resolution (A) gradient': [np.nan]*df_len,
                            'Resolution (A) r^2': [np.nan]*df_len,
                            'Resolution (A) spearman cc': [np.nan]*df_len,
                            'Rwork gradient': [np.nan]*df_len,
                            'Rwork r^2': [np.nan]*df_len,
                            'Rwork spearman cc': [np.nan]*df_len,
                            'Rfree gradient': [np.nan]*df_len,
                            'Rfree r^2': [np.nan]*df_len,
                            'Rfree spearman cc': [np.nan]*df_len,
                            'Size (Da) gradient': [np.nan]*df_len,
                            'Size (Da) r^2': [np.nan]*df_len,
                            'Size (Da) spearman cc': [np.nan]*df_len,
                            'Temperature (K) gradient': [np.nan]*df_len,
                            'Temperature (K) r^2': [np.nan]*df_len,
                            'Temperature (K) spearman cc': [np.nan]*df_len,
                            'Percent Glu and Asp gradient': [np.nan]*df_len,
                            'Percent Glu and Asp r^2': [np.nan]*df_len,
                            'Percent Glu and Asp spearman cc': [np.nan]*df_len,
                            'Num Glu and Asp gradient': [np.nan]*df_len,
                            'Num Glu and Asp r^2': [np.nan]*df_len,
                            'Num Glu and Asp spearman cc': [np.nan]*df_len,
                            'Num terminal O atoms gradient': [np.nan]*df_len,
                            'Num terminal O atoms r^2': [np.nan]*df_len,
                            'Num terminal O atoms spearman cc': [np.nan]*df_len,
                            'Relative RSRZ percentile gradient': [np.nan]*df_len,
                            'Relative RSRZ percentile r^2': [np.nan]*df_len,
                            'Relative RSRZ percentile spearman cc': [np.nan]*df_len})

for row in range(df_len):
    pdb = df['PDB code'][row]
    sub_df_stats['PDB code'][row] = pdb
    resolution = df['Resolution (A)'][row]
    resolution_array = np.array(df['Resolution (A)'].tolist())
    print(row, pdb, resolution)

    closest_structure_indices = []
    closest_structure_resolutions = []
    for num in range(1000):
        index = (np.abs(resolution_array-resolution)).argmin()
        val = resolution_array[index]
        closest_structure_indices.append(index)
        closest_structure_resolutions.append(val)
        resolution_array[index] = np.inf

    min_resolution = min(closest_structure_resolutions)
    max_resolution = max(closest_structure_resolutions)

    additional_structure_indices = []
    for index, num in np.ndenumerate(resolution_array):
        if float(num) == float(min_resolution) or float(num) == float(max_resolution):
            additional_structure_indices.append(index[0])

    structure_indices = closest_structure_indices + additional_structure_indices
    sub_df = df.iloc[structure_indices]
    sub_df = sub_df.reset_index(drop=True)
    sub_df.to_csv('Resolution_filter/{}_nearest_resolution_structures.csv'.format(pdb), index=False)

    for prop in ['Resolution (A)', 'Rwork', 'Rfree', 'Temperature (K)', 'Size (Da)', 'Percent Glu and Asp', 'Num Glu and Asp', 'Num terminal O atoms', 'Relative RSRZ percentile']:
        bnet_list = sub_df['Bnet'].tolist()
        prop_list = sub_df[prop].tolist()

        bnet_list = [bnet_list[index] for index, val in enumerate(prop_list) if not np.isnan(val)]
        prop_list = [prop_list[index] for index, val in enumerate(prop_list) if not np.isnan(val)]

        gradient, intercept, r2, p_val, std_error = stats.linregress(prop_list, bnet_list)
        spearman_cc = stats.spearmanr(prop_list, bnet_list)[0]
        sub_df_stats['{} gradient'.format(prop)][row] = gradient
        sub_df_stats['{} r^2'.format(prop)][row] = r2
        sub_df_stats['{} spearman cc'.format(prop)][row] = spearman_cc

stats_df = pd.DataFrame(sub_df_stats)
stats_df.to_csv('Resolution_window_PDB_file_properties_filtered_stats.csv', index=False)




# Absolute and relative RSRZ scores
df = pd.read_csv('PDB_file_properties.csv')

for row in range(df.shape[0]):
    pdb = str(df['PDB code'][row])
    if '+' in pdb:
        pdb = str(pdb.replace('.00E+', 'E'))
        pdb = str(pdb.replace('1.20E+09', '12E8'))

        pdb_report = requests.get('http://www.ebi.ac.uk/pdbe/entry-files/download/{}_validation.xml'.format(pdb.lower())).text
        try:
            abs_rsrz = float(pdb_report.split(' absolute-percentile-percent-RSRZ-outliers=')[1].split('"')[1])
        except:
            abs_rsrz = np.nan
        try:
            rel_rsrz = float(pdb_report.split(' relative-percentile-percent-RSRZ-outliers=')[1].split('"')[1])
        except:
            rel_rsrz = np.nan

        df.loc[row, 'PDB code'] = pdb
        df.loc[row, 'Absolute RSRZ percentile'] = abs_rsrz
        df.loc[row, 'Relative RSRZ percentile'] = rel_rsrz
