
import os
import requests
import shutil
import sys
import numpy as np
import pandas as pd

res_dict = {'ALA': 89.1,
            'ARG': 174.2,
            'ASN': 132.1,
            'ASP': 133.1,
            'CYS': 121.2,
            'GLN': 146.2,
            'GLU': 147.1,
            'GLY': 75.1,
            'HIS': 155.2,
            'ILE': 131.2,
            'LEU': 131.2,
            'LYS': 146.2,
            'MET': 149.2,
            'MSE': 196.1,
            'PHE': 165.2,
            'PRO': 115.1,
            'SER': 105.1,
            'THR': 119.1,
            'TRP': 204.2,
            'TYR': 181.2,
            'VAL': 117.1}

with open('../Protein_and_NA_PDB_IDs/Protein_PDB_IDs.txt', 'r') as f:
    protein_pdbs = f.readlines()
    protein_pdbs = [pdb.strip('\n') for pdb in protein_pdbs if len(pdb) == 5]

with open('../Protein_and_NA_PDB_IDs/NA_PDB_IDs.txt', 'r') as f:
    na_pdbs = f.readlines()
    na_pdbs = [pdb.strip('\n') for pdb in na_pdbs if len(pdb) == 5]

with open('../Protein_and_NA_PDB_IDs/XFEL_PDB_IDs.txt', 'r') as f:
    xfel_pdbs = f.readlines()
    xfel_pdbs = [pdb.strip('\n') for pdb in xfel_pdbs if len(pdb) == 5]

if bool(set(protein_pdbs) & set(na_pdbs)):
    sys.exit('Overlap between protein only and nucleic acid containing PDB'
             'accession code lists')

with open('../PDB_file_properties.csv', 'w') as f:
    f.write('PDB code' + ','
            'Resolution (A)' + ','
            'Rwork' + ','
            'Rfree' + ','
            'Temperature (K)' + ','
            'Size (kDa)' + ','
            'Contains NA?' + ','
            'Num Glu and Asp' + ','
            '% Glu and Asp' + ','
            'Electron density server response' + ','
            'Wilson plot B-factor (A^2)' + ','
            'Non-canonical aas count' + ','
            'PDBs no reprocessing (all conformers)' + ','
            'Automatically reprocessed PDBs (all conformers)' + ','
            'PDBs for manual reprocessing (all conformers)' + ','
            'PDBs no reprocessing (asp and glu conformers)' + ','
            'Automatically reprocessed PDBs (asp and glu conformers)' + ','
            'PDBs for manual reprocessing (asp and glu conformers)' + '\n')

with open('../PDBe_search_XFEL.csv', 'w') as f:
    f.write('PDB code' + ','
            'Resolution (A)' + ','
            'Rwork' + ','
            'Rfree' + ','
            'Temperature (K)' + ','
            'Size (kDa)' + ','
            'Contains NA?' + ','
            'Num Glu and Asp' + ','
            '% Glu and Asp' + ','
            'Electron density server response' + ','
            'Wilson plot B-factor (A^2)' + ','
            'Non-canonical aas count' + ','
            'PDBs no reprocessing (all conformers)' + ','
            'Automatically reprocessed PDBs (all conformers)' + ','
            'PDBs for manual reprocessing (all conformers)' + ','
            'PDBs no reprocessing (asp and glu conformers)' + ','
            'Automatically reprocessed PDBs (asp and glu conformers)' + ','
            'PDBs for manual reprocessing (asp and glu conformers)' + '\n')

if os.path.isdir('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/All_conformers/'):
    shutil.rmtree('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/All_conformers/')
os.makedirs('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/All_conformers/')

if os.path.isdir('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/Asp_Glu_conformers/'):
    shutil.rmtree('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/Asp_Glu_conformers/')
os.makedirs('/Volumes/Seagate_Backup_Plus_Drive_1TB/RABDAM/Asp_Glu_conformers/')

pdb_folder = '/Volumes/Seagate_Backup_Plus_Drive/pdb/asymmetric_unit/'
structure_count = 0
for middle_chars in os.listdir(pdb_folder):
    if len(middle_chars) == 2:
        for pdb_file_path in os.listdir('{}{}'.format(pdb_folder, middle_chars)):
            pdb_code = pdb_file_path.split('/')[-1][0:4].upper()

            structure_count += 1
            print('{}: {}'.format(pdb_code, structure_count))

            # Determines whether the structure contains nucleic acids
            contains_na = ''
            if pdb_code in protein_pdbs:
                contains_na = False
            elif pdb_code in na_pdbs:
                contains_na = True

            with open('{}{}/{}'.format(pdb_folder, middle_chars, pdb_file_path), 'r') as f:
                pdb_lines = f.readlines()

            # Initialises PDB file properties
            resolution = 0.0
            rwork = 0.0
            rfree = 0.0
            temperature = 0.0
            seqres_lines = []
            chain_res = []
            chain_id = ''
            xray = True
            electron_density = True
            wilson_b = 0.0

            ss_bonds = []
            all_new_pdb_lines = []
            asp_glu_new_pdb_lines = []
            all_reprocessed = False
            asp_glu_reprocessed = False

            all_raw_pdbs = ''
            all_reprocessed_pdbs = ''
            all_unprocessable_pdbs = ''
            asp_glu_raw_pdbs = ''
            asp_glu_reprocessed_pdbs = ''
            asp_glu_unprocessable_pdbs = ''

            # Checks that reliable electron density map can be downloaded from
            # the Uppsala Electron Density Server
            server_response = requests.get('http://eds.bmc.uu.se/cgi-bin/eds/'
                                           'uusfs?pdbCode={}'.format(pdb_code.lower())).text
            if 'sorry' in server_response:
                electron_density = False
            else:
                server_response_lines = server_response.split('\n')
                for index, line in enumerate(server_response_lines):
                    if 'Wilson plot B-factor' in line:
                        try:
                            wilson_b = float(server_response_lines[index+1].split('>')[1].split()[0])
                        except ValueError:
                            pass

            for index, line in enumerate(pdb_lines):
                line = line.strip('\n')
                no_whtspc_line = line.replace(' ', '')

                # Checks ATOM records for sub-1 occupancy non-alternate
                # conformers and increases their occcupancy to 1
                if (
                        line[0:6].strip() in ['ATOM', 'HETATM']
                    and line[17:20] in [res for chain in seqres_lines for res in chain]
                    and float(line[54:60]) != 1.0
                    and line[16:17] == ' '
                ):
                    all_reprocessed = True

                    all_new_pdb_line = line[:54] + '  1.00' + line[60:]
                    all_new_pdb_lines.append(all_new_pdb_line)
                else:
                    all_new_pdb_lines.append(line)

                # Checks ATOM records for non-alternate conformer sub-1
                # occupancy disulphide bonds and aspartate and glutamate
                # residues, and increases their occupancy to 1
                if (
                        line[0:4] == 'ATOM'
                    and line[17:20] in ['ASP', 'GLU']
                    and float(line[54:60]) != 1.0
                    and line[16:17] == ' '
                ):
                    asp_glu_reprocessed = True

                    asp_glu_new_pdb_line = line[:54] + '  1.00' + line[60:]
                    asp_glu_new_pdb_lines.append(asp_glu_new_pdb_line)
                elif (
                        line[0:4] == 'ATOM'
                    and line[17:20] == 'CYS'
                    and line[21:27].replace(' ', '') in ss_bonds
                    and float(line[54:60]) != 1.0
                    and line[16:17] == ' '
                ):
                    asp_glu_reprocessed = True

                    asp_glu_new_pdb_line = line[:54] + '  1.00' + line[60:]
                    asp_glu_new_pdb_lines.append(asp_glu_new_pdb_line)
                else:
                    asp_glu_new_pdb_lines.append(line)

                # Checks that structure is an X-ray crystal structure
                if (
                        no_whtspc_line.startswith('EXPDTA')
                    and not any(x in no_whtspc_line for x in
                            ['XRAYDIFFRACTION', 'X-RAYDIFFRACTION'])
                ):
                        xray = False
                        break

                # Finds structure resolution, and checks that it is < 5 Angstroms
                elif (
                        no_whtspc_line.startswith('REMARK2')
                    and 'RESOLUTION' in no_whtspc_line
                    and 'ANGSTROM' in no_whtspc_line
                ):
                    try:
                        resolution = float(line[23:30])
                        if resolution > 5:
                            resolution = 0.0
                            break
                    except ValueError:
                        resolution = 0.0
                        break

                # Finds structure Rwork
                elif (
                        no_whtspc_line.startswith('REMARK3RVALUE')
                    and any(x in no_whtspc_line for x in
                            ['(WORKINGSET):', '(WORKINGSET,NOCUTOFF):'])
                ):
                    sub_lines = no_whtspc_line.split(':')
                    try:
                        rwork = float(sub_lines[-1])
                    except ValueError:
                        rwork = 0.0
                        break

                # Finds structure Rfree
                elif (
                    any (x in no_whtspc_line for x in
                         ['REMARK3FREERVALUE:', 'REMARK3FREERVALUE(NOCUTOFF):'])
                ):
                    sub_lines = no_whtspc_line.split(':')
                    try:
                        rfree = float(sub_lines[-1])
                    except ValueError:
                        rfree = 0.0
                        break

                # Finds temperature of data collection
                elif no_whtspc_line.startswith('REMARK200TEMPERATURE(KELVIN):'):
                    sub_lines = no_whtspc_line.split(':')
                    try:
                        temperature = float(sub_lines[-1])
                    except ValueError:
                        temperature = 0.0
                        break

                # Finds residues in asymmetric unit
                elif line.startswith('SEQRES'):
                    res = line[19:].split()
                    new_chain_id = line[11:12]
                    if new_chain_id != chain_id:
                        if pdb_lines[index-1][0:6] == 'SEQRES':
                            seqres_lines.append(chain_res)
                            chain_res = []
                    chain_id = new_chain_id
                    chain_res += res

                elif line[0:6] != 'SEQRES' and pdb_lines[index-1][0:6] == 'SEQRES':
                    seqres_lines.append(chain_res)

                # Makes list of disulphide bonds
                elif line[0:6] == 'SSBOND':
                    ss_bonds.append(line[15:22].replace(' ', ''))
                    ss_bonds.append(line[29:36].replace(' ', ''))

            # Counts number of glutamate and aspartate residues and size of
            # structure
            count = 0
            size = 0.0
            glu_asp_count = 0.0
            non_canonical_aa_count = 0
            for chain in seqres_lines:
                sub_count = 0
                for res in chain:
                    # Excludes non-canonical amino acids (including nucleic acids)
                    if res in list(res_dict.keys()):
                        count += 1
                        sub_count += 1

                        # Finds structure size
                        size += res_dict[res]

                        # Finds structure Glu and Asp composition
                        if res in ['GLU', 'ASP']:
                            glu_asp_count += 1
                    else:
                        non_canonical_aa_count += 1

                size -= ((sub_count - 1) * 18.015)

            try:
                glu_asp_percent = ((glu_asp_count / count) * 100)
            except ZeroDivisionError:
                glu_asp_percent = 0.0
                pass

            # Determines whether structure is suitable for Bnet calculation or
            # requires reprocessing (all residues considered)
            if len(all_new_pdb_lines) > 0:
                all_atom_ids = []
                all_conformers = [[] for i in range(len(all_new_pdb_lines))]
                all_occupancy = [[] for i in range(len(all_new_pdb_lines))]

                # Lists alternate conformer identities and occupancies (per
                # atom)
                for line in all_new_pdb_lines:
                    if (
                            line[0:6].strip() in ['ATOM', 'HETATM']
                        and line[17:20] in [res for chain in seqres_lines for res in chain]
                        and float(line[54:60]) != 1.0
                    ):
                        atom_id = '{}_{}'.format(line[21:27].replace(' ', ''), line[12:16].strip())
                        if not atom_id in all_atom_ids:
                            all_atom_ids.append(atom_id)
                        index = all_atom_ids.index(atom_id)
                        if not line[16:17] in all_conformers[index]:
                            all_conformers[index].append(line[16:17])
                            all_occupancy[index].append(float(line[54:60]))

                # Checks if alternate conformer occupancies sum to 1
                all_unprocessable = False
                for index, atom_id in enumerate(all_atom_ids):
                    if atom_id.split('_')[0] in ss_bonds:
                        all_unprocessable = True

                    if np.sum(all_occupancy[index]) != 1.0:
                        all_unprocessable = True

                # Determines if PDB file needs reprocessing
                if all_reprocessed is False and all_unprocessable is False:
                    all_raw_pdbs = True
                elif all_reprocessed is True and all_unprocessable is False:
                    all_reprocessed_pdbs = True
                    with open(
                        '/Volumes/Seagate_Backup_Plus_Drive_1TB/sub_1_occ_PDB'
                        '_files_automatically_reprocessed/All_conformers/'
                        '{}_automatically_reprocessed.pdb'.format(pdb_code), 'w'
                    ) as f:
                        for line in all_new_pdb_lines:
                            f.write('{}\n'.format(line))
                else:
                    all_unprocessable_pdbs = True

            # Determines whether structure is suitable for Bnet calculation or
            # requires reprocessing (cysteine, aspartate and glutamate residues
            # considered only)
            if len(asp_glu_new_pdb_lines) > 0:
                asp_glu_atom_ids = []
                asp_glu_conformers = [[] for i in range(len(asp_glu_new_pdb_lines))]
                asp_glu_occupancy = [[] for i in range(len(asp_glu_new_pdb_lines))]

                # Lists alternate conformer identities and occupancies of
                # aspartate, glutamate and cysteine residues
                for line in asp_glu_new_pdb_lines:
                    if (
                        (
                        line[0:4] == 'ATOM'
                        and line[17:20] in ['ASP', 'GLU']
                        and float(line[54:60]) != 1.0
                        ) or (
                        line[0:4] == 'ATOM'
                        and line[17:20] == 'CYS'
                        and line[21:27].replace(' ', '') in ss_bonds
                        and float(line[54:60]) != 1.0
                        )
                    ):
                        atom_id = '{}_{}'.format(line[21:27].replace(' ', ''), line[12:16].strip())
                        if not atom_id in asp_glu_atom_ids:
                            asp_glu_atom_ids.append(atom_id)
                        index = asp_glu_atom_ids.index(atom_id)
                        if not line[16:17] in asp_glu_conformers[index]:
                            asp_glu_conformers[index].append(line[16:17])
                            asp_glu_occupancy[index].append(float(line[54:60]))

                # Checks if alternate conformer occupancies sum to 1
                asp_glu_unprocessable = False
                for index, atom_id in enumerate(asp_glu_atom_ids):
                    if atom_id.split('_')[0] in ss_bonds:
                        asp_glu_unprocessable = True

                    if np.sum(asp_glu_occupancy[index]) != 1.0:
                        asp_glu_unprocessable = True

                # Determines if PDB file needs reprocessing
                if asp_glu_reprocessed is False and asp_glu_unprocessable is False:
                    asp_glu_raw_pdbs = True
                elif asp_glu_reprocessed is True and asp_glu_unprocessable is False:
                    asp_glu_reprocessed_pdbs = True
                    with open(
                        '/Volumes/Seagate_Backup_Plus_Drive_1TB/sub_1_occ_PDB'
                        '_files_automatically_reprocessed/Asp_Glu_conformers/'
                        '{}_automatically_reprocessed.pdb'.format(pdb_code), 'w'
                    ) as f:
                        for line in asp_glu_new_pdb_lines:
                            f.write('{}\n'.format(line))
                else:
                    asp_glu_unprocessable_pdbs = True

            # Summarises PDB features in csv file
            if (
                   resolution == 0.0
                or rwork == 0.0
                or rfree == 0.0
                or temperature == 0.0
                or size == 0.0
                or contains_na == ''
                or glu_asp_count == 0.0
                or glu_asp_percent == 0.0
                or xray is False
            ):
                print('Unprocessed: {}'.format(pdb_code))
            else:
                with open('../PDB_file_properties.csv', 'a') as f:
                    f.write('{}'.format(pdb_code) + ','
                            '{}'.format(resolution) + ','
                            '{}'.format(rwork) + ','
                            '{}'.format(rfree) + ','
                            '{}'.format(temperature) + ','
                            '{}'.format(size) + ','
                            '{}'.format(contains_na) + ','
                            '{}'.format(glu_asp_count) + ','
                            '{}'.format(glu_asp_percent) + ','
                            '{}'.format(electron_density) + ','
                            '{}'.format(wilson_b) + ','
                            '{}'.format(non_canonical_aa_count) + ','
                            '{}'.format(all_raw_pdbs) + ','
                            '{}'.format(all_reprocessed_pdbs) + ','
                            '{}'.format(all_unprocessable_pdbs) + ','
                            '{}'.format(asp_glu_raw_pdbs) + ','
                            '{}'.format(asp_glu_reprocessed_pdbs) + ','
                            '{}'.format(asp_glu_unprocessable_pdbs) + '\n')

                if pdb_code in xfel_pdbs:
                    with open('../PDBe_search_XFEL.csv', 'a') as f:
                        f.write('{}'.format(pdb_code) + ','
                                '{}'.format(resolution) + ','
                                '{}'.format(rwork) + ','
                                '{}'.format(rfree) + ','
                                '{}'.format(temperature) + ','
                                '{}'.format(size) + ','
                                '{}'.format(contains_na) + ','
                                '{}'.format(glu_asp_count) + ','
                                '{}'.format(glu_asp_percent) + ','
                                '{}'.format(electron_density) + ','
                                '{}'.format(wilson_b) + ','
                                '{}'.format(non_canonical_aa_count) + ','
                                '{}'.format(all_raw_pdbs) + ','
                                '{}'.format(all_reprocessed_pdbs) + ','
                                '{}'.format(all_unprocessable_pdbs) + ','
                                '{}'.format(asp_glu_raw_pdbs) + ','
                                '{}'.format(asp_glu_reprocessed_pdbs) + ','
                                '{}'.format(asp_glu_unprocessable_pdbs) + '\n')

bnet_df = pd.read_csv('/Volumes/Seagate_Backup_Plus_Drive_1TB/Logfiles/Bnet_protein.csv')
all_stats_df = pd.read_csv('../PDB_file_properties.csv')
xfel_stats_df = pd.read_csv('../PDBe_search_XFEL.csv')
all_bnet_values = ['']*all_stats_df.shape[0]
xfel_bnet_values = ['']*xfel_stats_df.shape[0]

for bnet_index, pdb in enumerate(bnet_df['PDB'].tolist()):
    bnet = bnet_df['Bnet'][bnet_index]

    try:
        pdb_index = all_stats_df['PDB code'].tolist().index(pdb)
        all_bnet_values[pdb_index] = bnet
    except ValueError:
        pass

    try:
        pdb_index = xfel_stats_df['PDB code'].tolist().index(pdb)
        xfel_bnet_values[pdb_index] = bnet
    except ValueError:
        pass

all_bnet_values_df = pd.DataFrame({'Bnet': all_bnet_values})
all_df = pd.concat([all_stats_df, all_bnet_values_df], axis=1)
all_df.to_csv('../PDB_file_properties_plus_Bnet.csv')

xfel_bnet_values_df = pd.DataFrame({'Bnet': xfel_bnet_values})
xfel_df = pd.concat([xfel_stats_df, xfel_bnet_values_df], axis=1)
xfel_df.to_csv('../PDBe_search_XFEL_plus_Bnet.csv')
