
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
from Amino_acid_masses import amino_acid_masses
mass_dict = amino_acid_masses()


def check_rsrz_outliers(pdb_code):
    """
    Finds the Wilson B value, plus RSRZ outlier percentile values as compared
    to 1) all structures in the PDB and 2) structures of a similar resolution,
    from a structure's xml PDB validation report.
    *** Currently not in use to calculate RSRZ outlier values because RABDAM is
    being run on the PDB-REDO database, and these values as far as I can tell
    can only be obtained for the original PDB structures ***
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

    return wilson_b


def calc_bnet(pdb_code, cif_path, output_dir, rabdam_dir):
    """
    Runs RABDAM to calculate Bnet value
    """

    pdb_final = '{}_FINAL'.format(pdb_code.upper())
    bnet = np.nan

    os.system('echo "{},outputdir={},batchContinue=True,overwrite=True,PDT=7,'
              'windowSize=0.02,HETATM=Remove,proteinOrNucleicAcid=Protein"> '
              '{}/temp_rabdam_input.txt'.format(cif_path, output_dir, output_dir))
    os.system('python {}/rabdam/rabdam.py -i {}/temp_rabdam_input.txt -o '
              'bnet'.format(rabdam_dir, output_dir))
    os.remove('{}/temp_rabdam_input.txt'.format(output_dir))

    try:
        bnet_df = pd.read_pickle('{}/Logfiles/Bnet_protein.pkl'.format(output_dir))
        count = bnet_df['PDB'].tolist().count(pdb_final)
        if count != 1:
            pass
        else:
            try:
                index = bnet_df['PDB'].tolist().index(pdb_final)
                bnet = bnet_df['Bnet'][index]
            except ValueError:
                bnet = np.nan
                pass
    except FileNotFoundError:
        raise FileNotFoundError('{}/Logfiles/Bnet_protein.pkl doesn\'t exist'.format(output_dir))

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
        if '_pdbx_poly_seq_scheme.' in subsection:
            sub_lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
            prop_indices = {}
            prop_num = 0
            for line in sub_lines:
                if line.startswith('_pdbx_poly_seq_scheme'):
                    prop = line.split('.')[1].strip()
                    prop_indices[prop] = prop_num
                    prop_num += 1
                elif not line.startswith('_'):
                    line = line.split()
                    try:
                        chain_id = line[prop_indices['asym_id']]
                        res_id = line[prop_indices['mon_id']]
                        if not chain_id in list(seqres.keys()):
                            seqres[chain_id] = []
                        seqres[chain_id].append(res_id)
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
                            if (
                                    line[prop_indices['ptnr1_label_comp_id']] == 'CYS'
                                and line[prop_indices['ptnr2_label_comp_id']] == 'CYS'
                            ):
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
            # "null" = default weight of 1.0 used
            bfac_rest_str = line.split(':')[1].replace('\"', '').replace(',', '').strip()
            if bfac_rest_str.lower() == 'null':
                bfac_rest_weight = 1.0
            else:
                try:
                    bfac_rest_weight = float(bfac_rest_str)
                except ValueError:
                    pass
        elif '\"BREFTYPE\":' in line:
            bfac_model_raw = line.split(':')[1].replace('\"', '').replace(',', '').replace(' ', '')
            try:
                bfac_model = model_conv_dict[bfac_model_raw]
            except KeyError:
                pass
        # 0 = No TLS performed, otherwise denotes the number of TLS groups refined
        elif '\"NTLS\":' in line:
            tls_model_num = line.split(':')[1].replace('\"', '').replace(',', '')
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

    canonical_aas = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                     'THR', 'TRP', 'TYR', 'VAL']
    size = 0
    glu_asp_count = 0
    glu_asp_ld_count = 0
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
                if res in ['GLU', 'ASP', 'DGL', 'DAS']:
                    glu_asp_ld_count += 1
            if not res.upper() in canonical_aas:
                non_canonical_aa_count += 1

        size -= ((sub_count - 1) * 18.015)

    try:
        glu_asp_percent = ((glu_asp_count / count) * 100)
    except ZeroDivisionError:
        glu_asp_percent = 0

    # Converts to kDa
    size = size / 1000

    return (size, glu_asp_count, glu_asp_ld_count, glu_asp_percent,
            non_canonical_aa_count)


def find_peratom_properties(atom_lines, seqres, ss_bonds):
    """
    Finds structure properties in the PDB ATOM / HETATM records
    """

    sub_1_any = False
    sub_1_asp_glu = False
    sub_1_disulfide = False
    single_model = True

    seqres_list = [res for chain in seqres.values() for res in chain]
    ss_list = [res for bond in ss_bonds.values() for res in bond]

    alternate_conformers = {}
    terminal_o_atoms = []

    for line in atom_lines:
        # Checks ATOM records for sub-1 occupancy non-alternate conformers
        res_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                           line['pdbx_PDB_ins_code']])
        atom_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                            line['pdbx_PDB_ins_code'], line['label_comp_id'],
                            line['label_atom_id']])

        try:
            occupancy = float(line['occupancy'])
        except ValueError:
            sub_1_any = np.nan
            sub_1_asp_glu = np.nan
            sub_1_disulfide = np.nan
            alternate_conformers = {}
            break

        if (line['label_comp_id'] in seqres_list) and (occupancy != 1.0):
            if line['label_alt_id'] in ['.', '?']:
                sub_1_any = True
                if line['label_comp_id'] in ['ASP', 'GLU', 'DGL', 'DAS']:
                    sub_1_asp_glu = True
            else:
                if not atom_id in alternate_conformers:
                    alternate_conformers[atom_id] = {}
                alternate_conformers[atom_id][line['label_alt_id']] = occupancy

            # Checks for disulfide bonds with sub-1 occupancy
            if (res_id in ss_list) and (line['label_comp_id'] == 'CYS'):
                sub_1_disulfide = True

        # Finds terminal O atoms
        if (
                (line['label_comp_id'] in ['ASP', 'GLU', 'DGL', 'DAS'])
            and (line['label_atom_id'] in ['OD1', 'OD2', 'OE1', 'OE2'])
        ):
            terminal_o_atoms.append(atom_id)

    for line in atom_lines:
        # Checks only a single model present
        try:
            if float(line['pdbx_PDB_model_num']) != 1:
                single_model = False
        except ValueError:
            single_model = np.nan
            break

    # Checks that the occupancies of all alternate conformers sum to 1
    for atom_id, conformers_dict in alternate_conformers.items():
        occupancies = list(conformers_dict.values())
        if np.sum(occupancies) != 1.0:
            sub_1_any = True
            if atom_id.split('_')[3] in ['ASP', 'GLU', 'DGL', 'DAS']:
                sub_1_asp_glu = True

    # Calculates total number of terminal O atoms
    num_terminal_o_atoms = len(set(terminal_o_atoms))

    return (num_terminal_o_atoms, sub_1_any, sub_1_asp_glu,
            sub_1_disulfide, single_model)


def find_non_per_atom_b_factors(atom_lines, seqres):
    """
    Identifies structures with 20% or more residues with N, CA, C and O atoms
    with the same B-factor values
    """

    seqres_list = [res for chain in seqres.values() for res in chain]
    per_atom_b_factors = np.nan

    res_bfactors = {}
    for line in atom_lines:
        try:
            bfac = float(line['B_iso_or_equiv'])
        except ValueError:
            per_atom_b_factors = np.nan
            return per_atom_b_factors

        if line['label_comp_id'] in seqres_list:
            res_id = '_'.join([line['label_asym_id'], line['label_seq_id'],
                               line['pdbx_PDB_ins_code'], line['label_alt_id']])
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

    thresh = len(res_bfactors)*0.2
    if thresh < 3:
        thresh = 3
    if per_res_bfactor_count >= thresh:
        per_atom_b_factors = False
    else:
        per_atom_b_factors = True

    return per_atom_b_factors


def write_pdb_properties(num):
    """
    Writes properties of protein stuctures in the PDB circa 17th November 2020
    plus their Bnet values to a dataframe
    """

    print(num)

    cif_dir = '/home/shared/structural_bioinformatics/pdb_redo'
    work_dir = '/home/shared/structural_bioinformatics/RABDAM/PDB_parsing_output/{}/'.format(num)
    rabdam_dir = '/home/shared/structural_bioinformatics/RABDAM'

    with open(
        '{}/RABDAM_input_files/rabdam_protein_pdb_ids_part{}_201119.'
        'txt'.format(rabdam_dir, num), 'r'
    ) as f:
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
        all_pdbs_df = pd.read_pickle('{}/PDB_file_properties.pkl'.format(work_dir))
    else:
        all_pdbs_df = pd.DataFrame({'PDB code': [],
                                    'Resolution (A)': [],
                                    'Rwork': [],
                                    'Rfree': [],
                                    'Temperature (K)': [],
                                    'Size (kDa)': [],
                                    'Num Glu and Asp': [],
                                    'Num Glu and Asp (L and D)': [],
                                    '% Glu and Asp': [],
                                    'Num terminal O atoms': [],
                                    'Non-canonical aas count': [],
                                    'Wilson plot B-factor (A^2)': [],
                                    'Sub-1 occupancy (all conformers)': [],
                                    'Sub-1 occupancy (disulfide bonds)': [],
                                    'Sub-1 occupancy (asp and glu conformers)': [],
                                    'Single model': [],
                                    'Per-atom refined B-factors': [],
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
        cif_lines = None
        try:
            with open(
                '{}/{}/{}/{}_final.cif'.format(
                cif_dir, pdb_code[1:3].lower(), pdb_code.lower(), pdb_code.lower()
                ), 'r'
            ) as f:
                cif_lines = f.read().split('#')
        except FileNotFoundError:
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            continue

        # Opens PDB-REDO json file
        json_lines = None
        try:
            with open(
                '{}/{}/{}/data.json'.format(
                cif_dir, pdb_code[1:3].lower(), pdb_code.lower()
                ), 'r'
            ) as f:
                json_lines = f.read().split('\n')
        except FileNotFoundError:
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            continue

        # Generates dictionary of properties for each atom
        atom_lines = []
        for subsection in cif_lines:
            if '_atom_site.group_PDB' in subsection:
                lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
                prop_indices = {}
                prop_num = 0
                for line in lines:
                    if line.startswith('_atom_site.'):
                        prop = line.split('.')[1].strip()
                        prop_indices[prop] = prop_num
                        prop_num += 1
                    else:
                        if any(x in line for x in ['ATOM', 'HETATM']):
                            atom_props = {}
                            for prop, prop_index in prop_indices.items():
                                atom_props[prop] = line.split()[prop_index]
                            atom_lines.append(atom_props)

        # Finds structural properties of interest
        (resolution, rwork, rfree, temperature, seqres, ss_bonds
        ) = find_structurewide_properties(cif_lines)
        wilson_b = check_rsrz_outliers(pdb_code)
        (size, glu_asp_count, glu_asp_ld_count, glu_asp_percent,
         non_canonical_aa_count
        ) = find_aa_properties(seqres, mass_dict)
        (num_terminal_o_atoms, sub_1_any, sub_1_asp_glu, sub_1_disulfide,
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
            or np.isnan(single_model)
            or per_atom_b_factors is False
            or np.isnan(per_atom_b_factors)
        ):
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: {} unsuitable for Bnet calculation'.format(pdb_code))
            continue

        # Calculates Bnet
        cif_path = '{}/{}/{}/{}_final.cif'.format(
            cif_dir, pdb_code[1:3].lower(), pdb_code.lower(), pdb_code.lower()
        )
        bnet = calc_bnet(pdb_code, cif_path, work_dir, rabdam_dir)
        if np.isnan(bnet):
            with open('{}/Unprocessed_PDBs.txt'.format(work_dir), 'a') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to calculate Bnet for {}'.format(pdb_code))
            continue
        print('Bnet: {}'.format(bnet))

        # Summarises PDB features in dataframe
        print('Summarising info for {}'.format(pdb_code))
        indv_pdb_df = pd.DataFrame({'PDB code': [pdb_code],
                                    'Resolution (A)': [resolution],
                                    'Rwork': [rwork],
                                    'Rfree': [rfree],
                                    'Temperature (K)': [temperature],
                                    'Size (kDa)': [size],
                                    'Num Glu and Asp': [glu_asp_count],
                                    'Num Glu and Asp (L and D)': [glu_asp_ld_count],
                                    '% Glu and Asp': [glu_asp_percent],
                                    'Num terminal O atoms': [num_terminal_o_atoms],
                                    'Non-canonical aas count': [non_canonical_aa_count],
                                    'Wilson plot B-factor (A^2)': [wilson_b],
                                    'Sub-1 occupancy (all conformers)': [sub_1_any],
                                    'Sub-1 occupancy (disulfide bonds)': [sub_1_disulfide],
                                    'Sub-1 occupancy (asp and glu conformers)': [sub_1_asp_glu],
                                    'Single model': [single_model],
                                    'Per-atom refined B-factors': [per_atom_b_factors],
                                    'B-factor restraint weight': [bfac_rest_weight],
                                    'B-factor model': [bfac_model],
                                    'Number of TLS groups': [tls_model],
                                    'Bnet': [bnet]})
        all_pdbs_df = pd.concat([all_pdbs_df, indv_pdb_df], axis=0, ignore_index=True)
        all_pdbs_df.to_pickle('{}/PDB_file_properties.pkl'.format(work_dir))
        all_pdbs_df.to_excel('{}/PDB_file_properties.xlsx'.format(work_dir), index=False)


def calc_bnet_percentile():
    """
    Calculates bnet_percentile for PDB structure XXXX by comparing its Bnet
    value to the Bnet values of structures (min. 1000) of a similar resolution
    """

    work_dir = '/Users/ks17361/Lab_work_Elspeth_Garman/Papers/Bnet/PDB_analysis'
    # Filter to retain only those structures:
    # - Higher than 3.5 A resolution
    # - Temperature in the range of 80 - 120 K
    # - 20 or more Asp/Glu (L and D conformers) side-chain carnoxyl group oxygen atoms
    # - No disulfide bonds with multiple conformers
    # - No Asp/Glu residues with sub-1 occupancy (across all conformers)
    # - Single model
    # - Per-atom refined B-factors
    # - Non-flat B-factor model
    all_pdbs_df = pd.read_pickle('{}/PDB_file_properties_all_structures.pkl'.format(work_dir))
    filt_pdbs_df = all_pdbs_df[
          (all_pdbs_df['Resolution (A)'] <= 3.5)
        & (all_pdbs_df['Temperature (K)'] >= 80)
        & (all_pdbs_df['Temperature (K)'] <= 120)
        & (all_pdbs_df['Num terminal O atoms'] >= 20)
        & (all_pdbs_df['Sub-1 occupancy (disulfide bonds)'] == 0)
        & (all_pdbs_df['Sub-1 occupancy (asp and glu conformers)'] == 0)
        & (all_pdbs_df['Single model'] == 1)
        & (all_pdbs_df['Per-atom refined B-factors'] == 1)
        & (all_pdbs_df['B-factor model'].isin(['isotropic', 'anisotropic']))
    ]
    filt_pdbs_df = filt_pdbs_df.reset_index(drop=True)

    resolution_percentile = [np.nan]*filt_pdbs_df.shape[0]
    restraint_percentile = [np.nan]*filt_pdbs_df.shape[0]

    for row in range(filt_pdbs_df.shape[0]):
        pdb_code = filt_pdbs_df['PDB code'][row]
        resolution = filt_pdbs_df['Resolution (A)'][row]
        restraint = filt_pdbs_df['B-factor restraint weight'][row]
        bnet = filt_pdbs_df['Bnet'][row]
        if np.isnan(bnet) or np.isinf(bnet):
            raise ValueError(
                'WARNING: Not calculating Bnet percentile for {}'.format(pdb_code)
            )

        print('Calculating Bnet percentile for {} {}'.format(
            pdb_code, ((row+1) / filt_pdbs_df.shape[0])
        ))

        for percentile_tup in [
             [resolution, 'Resolution (A)', resolution_percentile],
             [restraint, 'B-factor restraint weight', restraint_percentile]
        ]:
            prop_val = percentile_tup[0]
            prop_name = percentile_tup[1]
            percentile = percentile_tup[2]

            array = copy.deepcopy(filt_pdbs_df[prop_name]).to_numpy()
            surr_struct_indices = []
            surr_struct_vals = []
            for num in range(1000):
                index = (np.abs(array-prop_val)).argmin()
                nearest_prop_val = array[index]
                surr_struct_indices.append(index)
                surr_struct_vals.append(nearest_prop_val)
                array[index] = np.inf

            min_val = min(surr_struct_vals)
            max_val = max(surr_struct_vals)
            for index, num in np.ndenumerate(array):
                if num == min_val or num == max_val:
                    surr_struct_indices.append(index[0])

            surr_struct_df = copy.deepcopy(filt_pdbs_df).iloc[surr_struct_indices].reset_index(drop=True)
            bnet_range = np.sort(surr_struct_df['Bnet'].to_numpy())
            bnet_percentile_list = np.where(bnet_range == bnet)
            if len(bnet_percentile_list) == 1:
                i = 0
                bnet_percentile = (bnet_percentile_list[i][0] + 1) / bnet_range.shape[0]
            elif len(bnet_percentile_list) > 1:
                if len(bnet_percentile_list) % 2 == 1:
                    i = (len(bnet_percentile_list) - 1) / 2
                    bnet_percentile = (bnet_percentile_list[i][0] + 1) / bnet_range.shape[0]
                elif len(bnet_percentile_list) % 2 == 0:
                    i = (len(bnet_percentile_list) - 2) / 2
                    j = len(bnet_percentile_list) / 2
                    average = (bnet_percentile_list[i][0] + bnet_percentile_list[j][0]) / 2
                    bnet_percentile = (average + 1) / bnet_range.shape[0]
            else:
                raise Exception(
                    'ERROR occured during calculation of Bnet percentile for '
                    '{}'.format(pdb_code)
                )
            percentile[row] = bnet_percentile

    bnet_percentile_df = pd.DataFrame({'Bnet percentile (resolution)': resolution_percentile,
                                       'Bnet percentile (restraint weight)': restraint_percentile})
    bnet_percentile_df = pd.concat(
        [filt_pdbs_df, bnet_percentile_df], axis=1
    ).reset_index(drop=True)
    bnet_percentile_df.to_pickle(
        '{}/PDB_file_properties_filtered_Bnet_percentile.pkl'.format(work_dir)
    )
    bnet_percentile_df.to_excel(
        '{}/PDB_file_properties_filtered_Bnet_percentile.xlsx'.format(work_dir), index=False
    )
