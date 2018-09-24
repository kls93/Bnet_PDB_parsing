
import pandas as pd
from collections import OrderedDict

pdb_codes_file = input('Specify csv file to which to add features: ')
bnet_df = pd.read_csv('PDB_file_properties.csv')
df = pd.read_csv(pdb_codes_file)

bnet_pdb_codes = bnet_df['PDB code'].tolist()
df_pdb_codes = [pdb for pdb in df['PDB ID'].tolist() if pdb != '']

resolution_values = ['']*df.shape[0]
rwork_values = ['']*df.shape[0]
rfree_values = ['']*df.shape[0]
temperature_values = ['']*df.shape[0]
size_values = ['']*df.shape[0]
glu_asp_percent_values = ['']*df.shape[0]
contains_na_values = ['']*df.shape[0]
bnet_values = ['']*df.shape[0]

for df_index, pdb in enumerate(df_pdb_codes):
    pdb = pdb.upper()
    try:
        bnet_index = bnet_pdb_codes.index(pdb)

        resolution = bnet_df['Resolution (A)'][bnet_index]
        rwork = bnet_df['Rwork'][bnet_index]
        rfree = bnet_df['Rfree'][bnet_index]
        temperature = bnet_df['Temperature (K)'][bnet_index]
        size = bnet_df['Size (kDa)'][bnet_index]
        glu_asp_percent = bnet_df['% Glu and Asp'][bnet_index]
        contains_na = bnet_df['Contains NA?'][bnet_index]
        bnet = bnet_df['Bnet'][bnet_index]

        resolution_values[df_index] = resolution
        rwork_values[df_index] = rwork
        rfree_values[df_index] = rfree
        temperature_values[df_index] = temperature
        size_values[df_index] = size
        glu_asp_percent_values[df_index] = glu_asp_percent
        contains_na_values[df_index] = contains_na
        bnet_values[df_index] = bnet

    except ValueError:
        pass

features_df = pd.DataFrame(OrderedDict({'Resolution (A)': resolution_values,
                                        'Rfree': rfree_values,
                                        'Rwork': rwork_values,
                                        'Temperature (K)': temperature_values,
                                        'Size (kDa)': size_values,
                                        '% Glu and Asp': glu_asp_percent_values,
                                        'Contains NA?': contains_na_values,
                                        'Bnet': bnet_values}))
df = pd.concat([df, features_df], axis=1)
df.to_csv('New_csv.csv')
