
import numpy as np
import pandas as pd
from collections import OrderedDict

file_pairs = {
    'Radiation_damage_series_Bnet_values.csv': '/Users/ks17361/Downloads/Rad_dam_series/',
    'PDBe_search_XFEL.csv': '/Users/ks17361/Downloads/XFEL/'
}

for orig_csv, reprocessed_bnet_csv in file_pairs.items():
    orig_df = pd.read_csv(orig_csv)
    new_bnet_df = pd.read_csv('{}Reprocessed/Logfiles/Bnet_protein.csv'.format(reprocessed_bnet_csv))

    with open('{}PDBs_no_reprocessing.txt'.format(reprocessed_bnet_csv), 'r') as f:
        raw_pdb_list = [pdb.strip('\n').upper() for pdb in f.readlines() if
                        len(pdb) == 5]
    with open('{}Automatically_reprocessed_PDBs.txt'.format(reprocessed_bnet_csv), 'r') as f:
        automatically_reprocessed_pdb_list = [pdb.strip('\n').upper() for pdb in
                                              f.readlines() if len(pdb) == 5]
    with open('{}PDBs_for_manual_reprocessing.txt'.format(reprocessed_bnet_csv), 'r') as f:
        manually_reprocessed_pdb_list = [pdb.strip('\n').upper() for pdb in
                                         f.readlines() if len(pdb) == 5]

    orig_pdbs = [pdb for pdb in orig_df['PDB ID'].tolist()]
    reprocessed_pdbs = [pdb for pdb in new_bnet_df['PDB'].tolist()]
    reprocessed_bnet = ['']*orig_df.shape[0]
    pdbs_no_reprocessing = ['']*orig_df.shape[0]
    automatically_reprocessed_pdbs = ['']*orig_df.shape[0]
    pdbs_for_manual_reprocessing = ['']*orig_df.shape[0]

    for orig_index, pdb in enumerate(orig_pdbs):
        if not pdb in ['', np.nan]:
            if pdb.upper() in raw_pdb_list:
                pdbs_no_reprocessing[orig_index] = 'Yes'
            if pdb.upper() in automatically_reprocessed_pdb_list:
                automatically_reprocessed_pdbs[orig_index] = 'Yes'
            if pdb.upper() in manually_reprocessed_pdb_list:
                pdbs_for_manual_reprocessing[orig_index] = 'Yes'

            try:
                bnet_index = reprocessed_pdbs.index('{}_REFMAC'.format(pdb.upper()))
                new_bnet = new_bnet_df['Bnet'][bnet_index]
                reprocessed_bnet[orig_index] = new_bnet
            except ValueError:
                pass

    reprocessed_bnet_df = pd.DataFrame(OrderedDict({'Re-refined Bnet': reprocessed_bnet,
                                                    'PDBs_no_reprocessing': pdbs_no_reprocessing,
                                                    'Automatically_reprocessed_PDBs': automatically_reprocessed_pdbs,
                                                    'PDBs_for_manual_reprocessing': pdbs_for_manual_reprocessing}))
    df = pd.concat([orig_df, reprocessed_bnet_df], axis=1)
    df.to_csv('{}_updated.csv'.format(orig_csv.split('.')[0]))
