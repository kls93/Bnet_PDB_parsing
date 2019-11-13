
# Functions to extract features and Bnet values for a list of PDB accession codes
# NOTE: No need to filter for nucleic acids because am only considering
# protein-only PDB codes, but would need to add in necessary filters if a
# different subset of PDB codes were used

import copy
import os
import requests
import shutil
import sys
import numpy as np
import pandas as pd


# Residue masses
mass_dict = {'ALA': 89.1,
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


def check_e_dens(pdb_code):
    """
    Checks that reliable electron density map can be downloaded from the
    Uppsala Electron Density Server
    """

    electron_density = True
    server_response = requests.get('http://eds.bmc.uu.se/cgi-bin/eds/'
                                   'uusfs?pdbCode={}'.format(pdb_code.lower())).text
    if 'sorry' in server_response:
        electron_density = False

    return electron_density


def check_rsrz_outliers(pdb_code):
    """
    Finds the Wilson B value, plus RSRZ outlier percentile values as compared
    to 1) all structures in the PDB and 2) structures of a similar resolution,
    from a structure's xml PDB validation report.
    """

    rsrz_outliers_all = np.nan
    rsrz_outliers_sim_res = np.nan
    wilson_b = np.nan

    pdb_info = requests.get(
        'http://files.rcsb.org/pub/pdb/validation_reports/{}/{}/{}_validation.xml.gz'.format(
        pdb_code.lower()[1:3], pdb_code.lower(), pdb_code.lower())
    )
    if pdb_info.status_code == 200:
        pdb_info = pdb_info.text
        try:
            rsrz_outliers_all = float(
                pdb_info.split('" absolute-percentile-percent-RSRZ-outliers=')[1].split('"')[1]
            )
        except (ValueError, IndexError):
            pass
        try:
            rsrz_outliers_sim_res = float(
                pdb_info.split('" relative-percentile-percent-RSRZ-outliers=')[1].split('"')[1]
            )
        except (ValueError, IndexError):
            pass
        try:
            wilson_b = float(pdb_info.split('" WilsonBestimate=')[1].split('"')[1])
        except (ValueError, IndexError):
            pass

    return rsrz_outliers_all, rsrz_outliers_sim_res, wilson_b


def calc_bnet(pdb_code):
    """
    Runs RABDAM to calculate Bnet value
    """

    bnet = np.nan

    os.system('echo "{},dir=/Volumes/Seagate_Backup_Plus_Drive/RABDAM_for_PDB/,'
              'batchContinue=True,overwrite=True"> temp_rabdam_input.txt'.format(pdb_code))
    os.system('python /Volumes/Seagate_Backup_Plus_Drive/RABDAM_for_PDB/RABDAM'
              '/rabdam/rabdam.py -i temp_rabdam_input.txt -o bnet')
    os.remove('temp_rabdam_input.txt')

    try:
        bnet_df = pd.read_pickle('/Volumes/Seagate_Backup_Plus_Drive/RABDAM_for_PDB'
                             '/Logfiles/Bnet_Protein.pkl')
        try:
            index = bnet_df['PDB'].tolist().index(pdb_code.upper())
            bnet = bnet_df['Bnet'][index]
            print(bnet)
        except ValueError:
            pass
    except FileNotFoundError:
        pass

    return bnet


def find_structurewide_properties(pdb_lines):
    """
    Finds structure properties in the PDB header
    """

    resolution = np.nan
    rwork = np.nan
    rfree = np.nan
    temperature = np.nan
    seqres = {}

    for subsection in pdb_lines:
        # Finds structure resolution
        if '_refine.ls_d_res_high' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_d_res_high'):
                    try:
                        resolution = float(line.replace('_refine.ls_d_res_high', ''))
                    except ValueError:
                        pass
                    break

        # Finds structure Rwork
        if '_refine.ls_R_factor_R_work' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_R_factor_R_work'):
                    try:
                        rwork = float(line.replace('_refine.ls_R_factor_R_work', ''))
                    except ValueError:
                        pass
                    break

        # Finds structure Rfree
        if '_refine.ls_R_factor_R_free' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_R_factor_R_free'):
                    try:
                        rfree = float(line.replace('_refine.ls_R_factor_R_free', ''))
                    except ValueError:
                        pass
                    break

        # Finds temperature of data collection
        if '_diffrn.ambient_temp' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_diffrn.ambient_temp'):
                    try:
                        temperature = float(line.replace('_diffrn.ambient_temp', ''))
                    except ValueError:
                        pass
                    break

        # Finds residues in asymmetric unit
        if '_entity_poly_seq.' in subsection:
            sub_lines = [line for line in subsection.split('\n') if not line.strip() == '']
            prop_indices = {}
            prop_num = 0
            for line in sub_lines:
                if line.startswith('_entity_poly_seq.'):
                    prop = line.split('.')[1].strip()
                    prop_indices[prop] = prop_num
                    prop_num += 1
                elif not any(line.startswith(x) for x in ['_', 'loop_']):
                    line = line.split()
                    i1 = prop_indices['entity_id']
                    if not line[i1] in seqres:
                        seqres[line[i1]] = []
                    i2 = prop_indices['mon_id']
                    seqres[line[i1]].append(line[i2].upper())

    print(resolution, rwork, rfree, temperature)

    return (resolution, rwork, rfree, temperature, seqres)


def find_aa_properties(seqres, mass_dict):
    """
    Calculates structure size and asp / glu content from seqres lines of PDB file
    """

    size = 0
    glu_asp_count = 0
    non_canonical_aa_count = 0

    # Counts number of glutamate and aspartate residues and size of structure
    # (excluding non-canonical amino acids other than MSE)
    count = 0
    for chain in seqres.keys():
        sub_count = 0
        for res in seqres[chain]:
            count += 1
            sub_count += 1

            if res in list(mass_dict.keys()):
                size += mass_dict[res]
                if res in ['GLU', 'ASP']:
                    glu_asp_count += 1
            else:
                non_canonical_aa_count += 1

        size -= ((sub_count - 1) * 18.015)

    try:
        glu_asp_percent = ((glu_asp_count / count) * 100)
    except ZeroDivisionError:
        glu_asp_percent = 0

    return size, glu_asp_count, glu_asp_percent, non_canonical_aa_count


def find_peratom_properties(atom_lines, seqres):
    """
    Finds structure properties in the PDB ATOM / HETATM records
    """

    all_reprocessed = False
    asp_glu_reprocessed = False
    single_model = True

    alternate_conformers = {}
    terminal_o_atoms = []

    for i, line in atom_lines.items():
        # Checks ATOM records for sub-1 occupancy non-alternate conformers
        if (
                line['label_comp_id'] in [res for chain in seqres.values() for res in chain]
            and float(line['occupancy']) != 1.0
        ):
            if line['label_alt_id'] == ' ':
                all_reprocessed = True
                if line['label_comp_id'] in ['ASP', 'GLU']:
                    asp_glu_reprocessed = True
            elif line['label_alt_id'] != ' ':
                atom_id = '_'.join([line['label_entity_id'], line['label_seq_id'],
                                    line['label_comp_id'], line['label_atom_id']])
                if not atom_id in alternate_conformers:
                    alternate_conformers[atom_id] = []
                alternate_conformers[atom_id].append(float(line['occupancy']))

        # Finds terminal O atoms
        if (
                line['label_comp_id'] in ['ASP', 'GLU']
            and line['label_atom_id'].strip() in ['OD1', 'OD2', 'OE1', 'OE2']
        ):
            atom_id = '_'.join([line['label_entity_id'], line['label_seq_id'],
                                line['label_comp_id'], line['label_atom_id']])
            terminal_o_atoms.append(atom_id)

        # Checks only a single model present
        if float(line['pdbx_PDB_model_num']) != 1:
            single_model = False

    # Checks that the occupancies of all alternate conformers sum to 1
    for atom_id in alternate_conformers.keys():
        if sum(alternate_conformers[atom_id]) != 1.0:
            all_reprocessed = True
            if atom_id.split('_')[2] in ['ASP', 'GLU']:
                asp_glu_reprocessed = True

    # Calculates total number of terminal O atoms
    num_terminal_o_atoms = len(set(terminal_o_atoms))

    return num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed, single_model


def find_non_per_atom_b_factors(atom_lines, seqres):
    """
    Identifies structures with 20% or more residues with N, CA, C and O atoms
    with the same B-factor values
    """

    res_bfactors = {}
    for i, line in atom_lines.items():
        if line['label_comp_id'] in [res for chain in seqres.values() for res in chain]:
            res_id = '_'.join([line['label_entity_id'], line['label_seq_id'],
                               line['label_alt_id'], line['label_comp_id']])
            atom_id = line['label_atom_id']
            if not res_id in res_bfactors:
                res_bfactors[res_id] = {}
            res_bfactors[res_id][atom_id] = float(line['B_iso_or_equiv'])

    per_res_bfactor_count = 0
    for res_id in res_bfactors:
        main_chain_bfactors = []

        for atom_id in res_bfactors[res_id]:
            if atom_id in ['N', 'CA', 'C', 'O']:
                main_chain_bfactors.append(res_bfactors[res_id][atom_id])

        if (
                (len(main_chain_bfactors) == 4)
            and (len(set(main_chain_bfactors)) == 1)
        ):
            per_res_bfactor_count += 1

    per_atom_b_factors = ''
    thresh = len(res_bfactors)*0.2
    if thresh < 3:
        thresh = 3
    if per_res_bfactor_count >= thresh:
        per_atom_b_factors = False
    else:
        per_atom_b_factors = True

    return per_atom_b_factors


def write_pdb_properties():
    """
    Writes properties of protein stuctures in the PDB circa 9th November 2019
    plus their Bnet values to a dataframe
    """

    with open('Protein_and_NA_PDB_IDs/Protein_PDB_IDs.txt', 'r') as f:
        protein_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    with open('Protein_and_NA_PDB_IDs/NA_PDB_IDs.txt', 'r') as f:
        na_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    with open('Protein_and_NA_PDB_IDs/XFEL_PDB_IDs.txt', 'r') as f:
        xfel_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    if bool(set(protein_pdbs) & set(na_pdbs)):
        sys.exit('Overlap between protein only and nucleic acid containing PDB'
                 'accession code lists')

    all_pdbs_df = pd.DataFrame({'PDB code': [],
                                'Resolution (A)': [],
                                'Rwork': [],
                                'Rfree': [],
                                'Temperature (K)': [],
                                'Size (kDa)': [],
                                'Num Glu and Asp': [],
                                '% Glu and Asp': [],
                                'Num terminal O atoms': [],
                                'Non-canonical aas count': [],
                                'Per-atom refined B-factors': [],
                                'Electron density server response': [],
                                'Wilson plot B-factor (A^2)': [],
                                'RSRZ relative to all structures': [],
                                'RSRZ relative to similar res structures': [],
                                'Reprocessing required? (all conformers)': [],
                                'Reprocessing required? (asp and glu conformers)': [],
                                'Bnet': []})

    if os.path.isdir('PDB_parsing_output'):
        shutil.rmtree('PDB_parsing_output')
    os.mkdir('PDB_parsing_output')
    with open('PDB_parsing_output/Progress_track.txt', 'w') as f:
        f.write('PDBs processed so far:\n')
    with open('PDB_parsing_output/Unprocessed_PDBs.txt', 'w') as f:
        f.write('PDBs unable to be processed:\n')

    structure_count = 0
    for pdb_code in protein_pdbs:
        # Marks progress
        structure_count += 1
        print('{}: {}'.format(pdb_code, (structure_count/len(protein_pdbs))))
        with open('PDB_parsing_output/Progress_track.txt', 'a') as f:
            f.write('{}\n'.format(pdb_code))

        # Downloads PDB file of asymmetric unit and extracts
        pdb_request = requests.get(
            'https://files.rcsb.org/download/{}.cif'.format(pdb_code.lower())
        )
        if pdb_request.status_code != 200:
            with open('PDB_parsing_output/Unprocessed_PDBs.txt', 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to process {}'.format(pdb_code))
            continue

        pdb_lines = pdb_request.text.split('#')
        atom_lines = {}
        for subsection in pdb_lines:
            if '_atom_site.group_PDB' in subsection:
                lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
                prop_indices = {}
                prop_num = 0
                for i1, line in enumerate(lines):
                    if line.startswith('_atom_site.'):
                        prop = line.split('.')[1].strip()
                        prop_indices[prop] = prop_num
                        prop_num += 1
                    elif line[0:6].strip() in ['ATOM', 'HETATM']:
                        atom_lines[i1] = {}
                        for prop in prop_indices.keys():
                            atom_lines[i1][prop] = line.split()[prop_indices[prop]]

        # Finds structure properties of interest
        electron_density = check_e_dens(pdb_code)
        (rsrz_outliers_all, rsrz_outliers_sim_res, wilson_b
        ) = check_rsrz_outliers(pdb_code)
        (resolution, rwork, rfree, temperature, seqres
        ) = find_structurewide_properties(pdb_lines)
        (size, glu_asp_count, glu_asp_percent, non_canonical_aa_count
        ) = find_aa_properties(seqres, mass_dict)
        (num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed, single_model
        ) = find_peratom_properties(atom_lines, seqres)
        per_atom_b_factors = find_non_per_atom_b_factors(atom_lines, seqres)
        if (
               np.isnan(resolution)
            or np.isnan(temperature)
            or size == 0
            or num_terminal_o_atoms == 0
            or single_model is False
            or per_atom_b_factors is False
        ):
            with open('PDB_parsing_output/Unprocessed_PDBs.txt', 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: {} unsuitable for Bnet calculation'.format(pdb_code))
            continue
        bnet = calc_bnet(pdb_code)
        if np.isnan(bnet):
            with open('PDB_parsing_output/Unprocessed_PDBs.txt', 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to calculate Bnet for {}'.format(pdb_code))
            continue

        # Summarises PDB features in csv file
        print('Summarising info for {}'.format(pdb_code))
        indv_pdb_df = pd.DataFrame({'PDB code': [pdb_code],
                                    'Resolution (A)': [resolution],
                                    'Rwork': [rwork],
                                    'Rfree': [rfree],
                                    'Temperature (K)': [temperature],
                                    'Size (kDa)': [size],
                                    'Num Glu and Asp': [glu_asp_count],
                                    '% Glu and Asp': [glu_asp_percent],
                                    'Num terminal O atoms': [num_terminal_o_atoms],
                                    'Non-canonical aas count': [non_canonical_aa_count],
                                    'Per-atom refined B-factors': [per_atom_b_factors],
                                    'Electron density server response': [electron_density],
                                    'Wilson plot B-factor (A^2)': [wilson_b],
                                    'RSRZ relative to all structures': [rsrz_outliers_all],
                                    'RSRZ relative to similar res structures': [rsrz_outliers_sim_res],
                                    'Reprocessing required? (all conformers)': [all_reprocessed],
                                    'Reprocessing required? (asp and glu conformers)': [asp_glu_reprocessed],
                                    'Bnet': [bnet]})
        all_pdbs_df = pd.concat([all_pdbs_df, indv_pdb_df], axis=0, ignore_index=True)
        all_pdbs_df.to_pickle('PDB_parsing_output/PDB_file_properties.pkl')

    return all_pdbs_df


def calc_bnet_percentile(all_pdbs_df):
    """
    Calculates bnet_percentile for PDB structure XXXX by comparing its Bnet
    value to the Bnet values of structures (min. 1000) of a similar resolution
    """

    bnet_percentile_df = [np.nan]*all_pdbs_df.shape[0]

    for row in range(all_pdbs_df.shape[0]):
        pdb_code = all_pdbs_df['PDB code'][row]
        resolution = all_pdbs_df['Resolution (A)'][row]
        bnet = all_pdbs_df['Bnet'][row]
        if np.isnan(bnet) or np.isinf(bnet):
            print('WARNING: Not calculating Bnet percentile for {}'.format(pdb_code))
            continue

        print('Calculating Bnet percentile for {} {}'.format(
            pdb_code, ((row+1) / all_pdbs_df.shape[0])
        ))

        resolution_array = copy.deepcopy(all_pdbs_df['Resolution (A)']).to_numpy()

        surr_struct_indices = []
        surr_struct_resns = []
        for num in range(1000):
            index = (np.abs(resolution_array-resolution)).argmin()
            nearest_resn = resolution_array[index]
            surr_struct_indices.append(index)
            surr_struct_resns.append(nearest_resn)
            resolution_array[index] = np.inf

        min_resolution = min(surr_struct_resns)
        max_resolution = max(surr_struct_resns)
        for index, num in np.ndenumerate(resolution_array):
            if num == min_resolution or num == max_resolution:
                surr_struct_indices.append(index[0])

        surr_struct_df = all_pdbs_df.iloc[surr_struct_indices].reset_index(drop=True)
        bnet_range = np.sort(surr_struct_df['Bnet'].to_numpy())
        bnet_percentile = (np.where(bnet_range == bnet)[0][0] + 1) / bnet_range.shape[0]
        bnet_percentile_df[row] = bnet_percentile

    bnet_percentile_df = pd.DataFrame({'Bnet percentile': bnet_percentile_df})
    bnet_percentile_df = pd.concat(
        [all_pdbs_df, bnet_percentile_df], axis=1
    )
    bnet_percentile_df.to_pickle('PDB_parsing_output/PDB_file_properties.pkl')
    bnet_percentile_df.to_csv('PDB_parsing_output/PDB_file_properties.csv', index=False)
    all_pdbs_df.to_pickle('PDB_parsing_output/Old_PDB_file_properties_just_in_case.pkl')

    return all_pdbs_df, bnet_percentile_df
