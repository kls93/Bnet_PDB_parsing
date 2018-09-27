
import numpy as np
import pandas as pd
from collections import OrderedDict

pdb_codes_file = input('Specify csv file to which to add features: ')
bnet_df = pd.read_csv('PDB_file_properties.csv')
df = pd.read_csv(pdb_codes_file)

bnet_pdb_codes = bnet_df['PDB code'].tolist()
df_pdb_codes = [pdb for pdb in df['PDB code'].tolist()]

resolution_values = ['']*df.shape[0]
rwork_values = ['']*df.shape[0]
rfree_values = ['']*df.shape[0]
temperature_values = ['']*df.shape[0]
size_values = ['']*df.shape[0]
contains_na_values = ['']*df.shape[0]
num_glu_asp_values = ['']*df.shape[0]
glu_asp_percent_values = ['']*df.shape[0]
e_dens_response_values = ['']*df.shape[0]
wilson_b_values = ['']*df.shape[0]
non_canonical_aa_count_values = ['']*df.shape[0]
all_no_reprocess_values = ['']*df.shape[0]
all_automatic_reprocess_values = ['']*df.shape[0]
all_manual_reprocess_values = ['']*df.shape[0]
glu_asp_no_reprocess_values = ['']*df.shape[0]
glu_asp_automatic_reprocess_values = ['']*df.shape[0]
glu_asp_manual_reprocess_values = ['']*df.shape[0]

# bnet_values = ['']*df.shape[0]

for df_index, pdb in enumerate(df_pdb_codes):
    if not pdb in ['', np.nan]:
        pdb = pdb.upper()
        try:
            bnet_index = bnet_pdb_codes.index(pdb)

            resolution = bnet_df['Resolution (A)'][bnet_index]
            rwork = bnet_df['Rwork'][bnet_index]
            rfree = bnet_df['Rfree'][bnet_index]
            temperature = bnet_df['Temperature (K)'][bnet_index]
            size = bnet_df['Size (kDa)'][bnet_index]
            contains_na = bnet_df['Contains NA?'][bnet_index]
            num_glu_asp = bnet_df['Num Glu and Asp'][bnet_index]
            glu_asp_percent = bnet_df['% Glu and Asp'][bnet_index]
            e_dens_response = bnet_df['Electron density server response'][bnet_index]
            wilson_b = bnet_df['Wilson plot B-factor (A^2)'][bnet_index]
            non_canonical_aa_count = bnet_df['Non-canonical aas count'][bnet_index]
            all_no_reprocess = bnet_df['PDBs no reprocessing (all conformers)'][bnet_index]
            all_automatic_reprocess = bnet_df['Automatically reprocessed PDBs (all conformers)'][bnet_index]
            all_manual_reprocess = bnet_df['PDBs for manual reprocessing (all conformers)'][bnet_index]
            glu_asp_no_reprocess = bnet_df['PDBs no reprocessing (asp and glu conformers)'][bnet_index]
            glu_asp_automatic_reprocess = bnet_df['Automatically reprocessed PDBs (asp and glu conformers)'][bnet_index]
            glu_asp_manual_reprocess = bnet_df['PDBs for manual reprocessing (asp and glu conformers)'][bnet_index]

            # bnet = bnet_df['Bnet'][bnet_index]

            resolution_values[df_index] = resolution
            rwork_values[df_index] = rwork
            rfree_values[df_index] = rfree
            temperature_values[df_index] = temperature
            size_values[df_index] = size
            contains_na_values[df_index] = contains_na
            num_glu_asp_values[df_index] = num_glu_asp
            glu_asp_percent_values[df_index] = glu_asp_percent
            e_dens_response_values[df_index] = e_dens_response
            wilson_b_values[df_index] = wilson_b
            non_canonical_aa_count_values[df_index] = non_canonical_aa_count
            all_no_reprocess_values[df_index] = all_no_reprocess
            all_automatic_reprocess_values[df_index] = all_automatic_reprocess
            all_manual_reprocess_values[df_index] = all_manual_reprocess
            glu_asp_no_reprocess_values[df_index] = glu_asp_no_reprocess
            glu_asp_automatic_reprocess_values[df_index] = glu_asp_automatic_reprocess
            glu_asp_manual_reprocess_values[df_index] = glu_asp_manual_reprocess

            # bnet_values[df_index] = bnet

        except ValueError:
            pass

features_df = pd.DataFrame(OrderedDict({'Resolution (A)': resolution_values,
                                        'Rfree': rfree_values,
                                        'Rwork': rwork_values,
                                        'Temperature (K)': temperature_values,
                                        'Size (kDa)': size_values,
                                        'Contains NA?': contains_na_values,
                                        'Num Glu and Asp': num_glu_asp_values,
                                        '% Glu and Asp': glu_asp_percent_values,
                                        'Electron density server response': e_dens_response_values,
                                        'Wilson plot B-factor (A^2)': wilson_b_values,
                                        'Non-canonical aas count': non_canonical_aa_count_values,
                                        'PDBs no reprocessing (all conformers)': all_no_reprocess_values,
                                        'Automatically reprocessed PDBs (all conformers)': all_automatic_reprocess_values,
                                        'PDBs for manual reprocessing (all conformers)': all_manual_reprocess_values,
                                        'PDBs no reprocessing (asp and glu conformers)': glu_asp_no_reprocess_values,
                                        'Automatically reprocessed PDBs (asp and glu conformers)': glu_asp_automatic_reprocess_values,
                                        'PDBs for manual reprocessing (asp and glu conformers)': glu_asp_manual_reprocess_values}))
df = pd.concat([df, features_df], axis=1)
df.to_csv('New_csv.csv')
