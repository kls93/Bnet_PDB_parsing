
import os
import pandas as pd

pdb_folder = '/Volumes/Seagate_Backup_Plus_Drive/pdb/asymmetric_unit/'
structure_count = 0
per_atom_refined_structures = {}

side_chain_atom_numbers = {'ALA': 1,
                           'ARG': 7,
                           'ASN': 4,
                           'ASP': 4,
                           'CYS': 2,
                           'GLN': 5,
                           'GLU': 5,
                           'HIS': 6,
                           'ILE': 4,
                           'LEU': 4,
                           'LYS': 5,
                           'MET': 4,
                           'PHE': 7,
                           'PRO': 3,
                           'SER': 2,
                           'THR': 3,
                           'TRP': 10,
                           'TYR': 8,
                           'VAL': 3}

for middle_chars in os.listdir(pdb_folder):
    if len(middle_chars) == 2:
        for pdb_file_path in os.listdir('{}{}'.format(pdb_folder, middle_chars)):
            pdb_code = pdb_file_path.split('/')[-1][0:4].upper()

            structure_count += 1
            print('{}: {}'.format(pdb_code, structure_count))

            per_atom_b_factors = False

            with open('{}{}/{}'.format(pdb_folder, middle_chars, pdb_file_path), 'r') as f:
                pdb_file_lines = f.readlines()

            res_bfactors = {}
            for line in pdb_file_lines:
                if line[0:6].strip() in ['ATOM']:
                    res_id = line[21:27].replace(' ', '') + '_' + line[16:20]
                    atom_id = line[12:16].strip()
                    if not res_id in res_bfactors:
                        res_bfactors[res_id] = {}
                    res_bfactors[res_id][atom_id] = float(line[60:66])

            for res_id in res_bfactors:
                main_chain_bfactors = []
                side_chain_bfactors = []

                res_name = res_id.split('_')[-1][1:]
                if res_name in side_chain_atom_numbers:
                    for atom_id in res_bfactors[res_id]:
                        if atom_id in ['N', 'CA', 'C', 'O']:
                            main_chain_bfactors.append(res_bfactors[res_id][atom_id])
                        else:
                            side_chain_bfactors.append(res_bfactors[res_id][atom_id])

                    if (
                        (    (len(main_chain_bfactors) == 4)
                         and (len(set(main_chain_bfactors)) == 1))
                        or
                        (    (len(side_chain_bfactors) == side_chain_numbers[res_name])
                         and (len(set(side_chain_bfactors)) == 1))
                    ):
                        per_atom_b_factors = True

            per_atom_refined_structures[pdb_code] = per_atom_b_factors

df = pd.DataFrame({'PDB code': list(per_atom_refined_structures.keys()),
                   'Per-atom refined?': list(per_atom_refined_structures.values())})
df.to_csv('../Per-atom_refined_PDBs.csv')
