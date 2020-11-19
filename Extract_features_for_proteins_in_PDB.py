
# Functions to extract features and Bnet values for a list of PDB accession codes
# NOTE: No need to filter for nucleic acids because am only considering
# protein-only PDB codes, but would need to add in necessary filters if a
# different subset of PDB codes were used

import copy
import os
import requests
import shutil
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


def calc_bnet(pdb_code, cif_path, output_dir, rabdam_dir):
    """
    Runs RABDAM to calculate Bnet value
    """

    bnet = np.nan

    os.system('echo "{},dir={},batchContinue=False,'
              'overwrite=True"> temp_rabdam_input.txt'.format(cif_path, output_dir))
    os.system('python {}/rabdam/rabdam.py -i temp_rabdam_input.txt -o bnet'.format(rabdam_dir))
    os.remove('temp_rabdam_input.txt')

    try:
        bnet_df = pd.read_pickle('{}/Logfiles/Bnet_Protein.pkl'.format(output_dir))
        try:
            index = bnet_df['PDB'].tolist().index(pdb_code.upper())
            bnet = bnet_df['Bnet'][index]
        except ValueError:
            pass
    except FileNotFoundError:
        raise FileNotFoundError('{}/Logfiles/Bnet_Protein.pkl doesn\'t exist'.format(output_dir))

    return bnet


def find_structurewide_properties(cif_lines):
    """
    Finds structure properties in the mmCIF file header
    """

    resolution = np.nan
    rwork = np.nan
    rfree = np.nan
    temperature = np.nan
    seqres = {}
    ss_bonds = {}

    for subsection in cif_lines:
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
            sub_lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
            prop_indices = {}
            prop_num = 0
            for line in sub_lines:
                if line.startswith('_entity_poly_seq.'):
                    prop = line.split('.')[1].strip()
                    prop_indices[prop] = prop_num
                    prop_num += 1
                elif not line.startswith('_'):
                    line = line.split()
                    try:
                        i1 = prop_indices['entity_id']
                        if not line[i1] in seqres:
                            seqres[line[i1]] = []
                        i2 = prop_indices['mon_id']
                        seqres[line[i1]].append(line[i2].upper())
                    except (KeyError, IndexError):
                        seqres = {}
                        break

        # Finds disulfide bonds
        if '_struct_conn.' in subsection:
            sub_lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
            prop_indices = {}
            prop_num = 0
            for line in sub_lines:
                if line.startswith('_struct_conn.'):
                    prop = line.split('.')[1].strip()
                    prop_indices[prop] = prop_num
                    prop_num += 1
                elif line.startswith('disulf'):
                    line = line.split()
                    try:
                        if line[prop_indices['conn_type_id']] == 'disulf':
                            bond_num = line[prop_indices['id']].replace('disulf', '')
                            res1 = (  line[prop_indices['ptnr1_label_asym_id']] + '_'
                                    + line[prop_indices['ptnr1_label_seq_id']] + '_'
                                    + line[prop_indices['pdbx_ptnr1_PDB_ins_code']])
                            res2 = (  line[prop_indices['ptnr2_label_asym_id']] + '_'
                                    + line[prop_indices['ptnr2_label_seq_id']] + '_'
                                    + line[prop_indices['pdbx_ptnr2_PDB_ins_code']])
                            ss_bonds[bond_num] = [res1, res2]
                    except (KeyError, IndexError):
                        ss_bonds = {}
                        break

    return (resolution, rwork, rfree, temperature, seqres, ss_bonds)

def find_b_factor_restraints(pdb_code, json_lines):
    """
    Finds the strength of the B-factor restraints applied in the final round of
    refinement of PDB-REDO structures
    """

    bfac_rest_weight = np.nan
    bfac_model = np.nan
    tls_model = np.nan

    model_conv_dict = {'ISOT': 'isotropic',
                       'ANISOT': 'anisotropic',
                       'OVER': 'flat'}

    for line in json_lines:
        if '\"BBEST\":' in line:
            try:
                bfac_rest_weight = float(line.split(':')[1].replace('\"', ''))
            except ValueError:
                pass
        elif '\"BREFTYPE\":' in line:
            bfac_model_raw = line.split(':')[1].replace('\"', '').replace(' ', '')
            try:
                bfac_model = model_conv_dict[bfac_model_raw]
            except KeyError:
                pass
        # 0 = No TLS performed, otherwise denotes the number of TLS groups refined
        elif '\"NTLS\":' in line:
            tls_model_num = line.split(':')[1]
            try:
                if int(tls_model_num) == float(tls_model_num):
                    tls_model = int(tls_model_num)
                else:
                    pass
            except ValueError:
                pass

    return (bfac_rest_weight, bfac_model, tls_model)


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

            if res.upper() in list(mass_dict.keys()):
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


def find_peratom_properties(atom_lines, seqres, ss_bonds):
    """
    Finds structure properties in the PDB ATOM / HETATM records
    """

    all_reprocessed = False
    asp_glu_reprocessed = False
    cys_reprocessed = False
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
                atom_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                                    line['pdbx_PDB_ins_code'], line['label_comp_id'],
                                    line['label_atom_id']])
                if not atom_id in alternate_conformers:
                    alternate_conformers[atom_id] = []
                alternate_conformers[atom_id].append(float(line['occupancy']))

            # Checks for disulfide bonds with sub-1 occupancy
            if (
                '_'.join([line['label_asym_id'], line['label_seq_id'], line['pdbx_PDB_ins_code']])
                in [res for bond in ss_bonds.values() for res in bond]
            ):
                cys_reprocessed = True

        # Finds terminal O atoms
        if (
                line['label_comp_id'] in ['ASP', 'GLU']
            and line['label_atom_id'] in ['OD1', 'OD2', 'OE1', 'OE2']
        ):
            atom_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                                line['pdbx_PDB_ins_code'], line['label_comp_id'],
                                line['label_atom_id']])
            terminal_o_atoms.append(atom_id)

        # Checks only a single model present
        if float(line['pdbx_PDB_model_num']) != 1:
            single_model = False

    # Checks that the occupancies of all alternate conformers sum to 1
    for atom_id in alternate_conformers.keys():
        if sum(alternate_conformers[atom_id]) != 1.0:
            all_reprocessed = True
            if atom_id.split('_')[3] in ['ASP', 'GLU']:
                asp_glu_reprocessed = True

    # Calculates total number of terminal O atoms
    num_terminal_o_atoms = len(set(terminal_o_atoms))

    return (num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed,
            cys_reprocessed, single_model)


def find_non_per_atom_b_factors(atom_lines, seqres):
    """
    Identifies structures with 20% or more residues with N, CA, C and O atoms
    with the same B-factor values
    """

    res_bfactors = {}
    for i, line in atom_lines.items():
        if line['label_comp_id'] in [res for chain in seqres.values() for res in chain]:
            res_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                               line['pdbx_PDB_ins_code'], line['label_alt_id'],
                               line['label_comp_id']])
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
    Writes properties of protein stuctures in the PDB circa 17th November 2020
    plus their Bnet values to a dataframe
    """

    cif_dir = '/home/shared/structural_bioinformatics/pdb_redo'
    work_dir = '/home/shared/structural_bioinformatics/RABDAM/PDB_parsing_output'
    rabdam_dir = 'home/shared/structural_bioinformatics/RABDAM'

    with open('{}/rabdam_protein_pdb_ids_201119.txt'.format(rabdam_dir), 'r') as f:
        protein_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    if not os.path.isfile('{}/Progress_track.txt'.format(work_dir)):
        with open('{}/Progress_track.txt'.format(work_dir), 'w') as f:
            f.write('PDBs processed so far:\n')
    if not os.path.isfile('{}/Unprocessed_PDBs.txt'.format(work_dir)):
        with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'w') as f:
            f.write('PDBs unable to be processed:\n')
    if os.path.isfile('{}/PDB_file_properties.pkl'.format(work_dir)):
        all_pdbs_df = pd.read_pickle('PDB_parsing_output/PDB_file_properties.pkl')
    else:
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
                                    'Sub-1 occupancy (all conformers)': [],
                                    'Sub-1 occupancy (disulfide bonds)': [],
                                    'Sub-1 occupancy (asp and glu conformers)': [],
                                    'B-factor restraint weight': [],
                                    'B-factor model': [],
                                    'Number of TLS groups': [],
                                    'Bnet': []})

    structure_count = 0
    for pdb_code in protein_pdbs:
        # Tracks progress
        structure_count += 1
        print('{}: {}'.format(pdb_code, (structure_count/len(protein_pdbs))))
        with open('{}/Progress_track.txt'.format(work_dir), 'a') as f:
            f.write('{}\n'.format(pdb_code))

        # Opens mmCIF PDB-REDO file
        try:
            with open(
                '{}/{}/{}/{}_final.cif'.format(cif_dir, pdb_code[1:3], pdb_code, pdb_code), 'r'
            ) as f:
                cif_lines = f.read().split('#')
        except FileNotFoundError:
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            continue

        # Opens PDB-REDO json file
        try:
            with open(
                '{}/{}/{}/data.json'.format(cif_dir, pdb_code[1:3], pdb_code), 'r'
            ) as f:
                json_lines = f.read().split('\n')
        except FileNotFoundError:
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            continue

        # Generates dictionary of properties for each atom
        atom_lines = {}
        for subsection in cif_lines:
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

        # Finds structural properties of interest
        electron_density = check_e_dens(pdb_code)
        (rsrz_outliers_all, rsrz_outliers_sim_res, wilson_b
        ) = check_rsrz_outliers(pdb_code)
        (resolution, rwork, rfree, temperature, seqres, ss_bonds
        ) = find_structurewide_properties(cif_lines)
        (size, glu_asp_count, glu_asp_percent, non_canonical_aa_count
        ) = find_aa_properties(seqres, mass_dict)
        (num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed, cys_reprocessed,
         single_model) = find_peratom_properties(atom_lines, seqres, ss_bonds)
        per_atom_b_factors = find_non_per_atom_b_factors(atom_lines, seqres)
        bfac_rest_weight, bfac_model, tls_model = find_b_factor_restraints(
            pdb_code, json_lines
        )

        # Filters out PDBs with unsuitable properties for Bnet calculation
        if (
               np.isnan(resolution)
            or np.isnan(temperature)
            or size == 0
            or num_terminal_o_atoms == 0
            or single_model is False
            or per_atom_b_factors is False
        ):
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: {} unsuitable for Bnet calculation'.format(pdb_code))
            continue

        # Calculates Bnet
        cif_path = '{}/{}/{}/{}_final.cif'.format(cif_dir, pdb_code[1:3], pdb_code, pdb_code)
        bnet = calc_bnet(pdb_code, cif_path, work_dir, rabdam_dir)
        if np.isnan(bnet):
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to calculate Bnet for {}'.format(pdb_code))
            continue

        # Summarises PDB features in dataframe
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
                                    'Sub-1 occupancy (all conformers)': [all_reprocessed],
                                    'Sub-1 occupancy (disulfide bonds)': [cys_reprocessed],
                                    'Sub-1 occupancy (asp and glu conformers)': [asp_glu_reprocessed],
                                    'B-factor restraint weight': [bfac_rest_weight],
                                    'B-factor model': [bfac_model],
                                    'Number of TLS groups': [tls_model],
                                    'Bnet': [bnet]})
        all_pdbs_df = pd.concat([all_pdbs_df, indv_pdb_df], axis=0, ignore_index=True)
        all_pdbs_df.to_pickle('PDB_parsing_output/PDB_file_properties.pkl')

    return all_pdbs_df


def calc_bnet_percentile(all_pdbs_df):
    """
    Calculates bnet_percentile for PDB structure XXXX by comparing its Bnet
    value to the Bnet values of structures (min. 1000) of a similar resolution
    """

    work_dir = '/home/shared/structural_bioinformatics/RABDAM/PDB_parsing_output'
    if not os.path.isfile('{}/Bnet_percentile_unprocessed_pdbs.txt'.format(work_dir)):
        with open('{}/Bnet_percentile_unprocessed_pdbs.txt'.format(work_dir), 'w') as f:
            f.write('Structures with a null / infinite Bnet value (which is '
                    'hence unsuitable for percentile calculations)\n')

    resolution_percentile = [np.nan]*all_pdbs_df.shape[0]
    restraint_percentile = [np.nan]*all_pdbs_df.shape[0]

    for row in range(all_pdbs_df.shape[0]):
        pdb_code = all_pdbs_df['PDB code'][row]
        resolution = all_pdbs_df['Resolution (A)'][row]
        restraint = all_pdbs_df['B-factor restraint weight'][row]
        bnet = all_pdbs_df['Bnet'][row]
        if np.isnan(bnet) or np.isinf(bnet):
            print('WARNING: Not calculating Bnet percentile for {}'.format(pdb_code))
            with open('{}/Bnet_percentile_unprocessed_pdbs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            continue

        print('Calculating Bnet percentile for {} {}'.format(
            pdb_code, ((row+1) / all_pdbs_df.shape[0])
        ))

        for percentile_tup in [
             [resolution, 'Resolution (A)', resolution_percentile],
             [restraint, 'B-factor restraint weight', restraint_percentile]
        ]:
            prop_num = percentile_tup[0]
            prop_name = percentile_tup[1]
            percentile = percentile_tup[2]

            array = copy.deepcopy(all_pdbs_df[prop_name]).to_numpy()
            surr_struct_indices = []
            surr_struct_vals = []
            for num in range(1000):
                index = (np.abs(array-prop_num)).argmin()
                nearest_prop_val = array[index]
                surr_struct_indices.append(index)
                surr_struct_vals.append(nearest_prop_val)
                array[index] = np.inf

            min_val = min(surr_struct_vals)
            max_val = max(surr_struct_vals)
            for index, num in np.ndenumerate(array):
                if num == min_val or num == max_val:
                    surr_struct_indices.append(index[0])

            surr_struct_df = copy.deepcopy(all_pdbs_df).iloc[surr_struct_indices].reset_index(drop=True)
            bnet_range = np.sort(surr_struct_df['Bnet'].to_numpy())
            bnet_percentile = (np.where(bnet_range == bnet)[0][0] + 1) / bnet_range.shape[0]
            percentile[row] = bnet_percentile

    bnet_percentile_df = pd.DataFrame({'Bnet percentile (resolution)': resolution_percentile,
                                       'Bnet percentile (restraint weight)': restraint_percentile})
    bnet_percentile_df = pd.concat(
        [all_pdbs_df, bnet_percentile_df], axis=1
    ).reset_index(drop=True)
    bnet_percentile_df.to_pickle('{}/PDB_file_properties_percentile.pkl'.format(work_dir))
    bnet_percentile_df.to_csv('{}/PDB_file_properties_percentile.csv'.format(work_dir), index=False)

    return bnet_percentile_df
