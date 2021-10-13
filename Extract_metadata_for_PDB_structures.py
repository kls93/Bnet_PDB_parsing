
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

def gen_pdb_list():
    pdb_list = ['2bhx',
                '2bi1',
                '2bi2',
                '2bi3',
                '2bi5',
                '2bi9',
                '2bia',
                '2j5k',
                '2j5q',
                '2j5r',
                '2yae',
                '2yaf',
                '2yah',
                '2yam',
                '2yao',
                '2yap',
                '2yaq',
                '2yar',
                '2ybh',
                '2ybi',
                '2ybj',
                '2ybl',
                '2ybm',
                '2ybn',
                '3ddz',
                '3de0',
                '3de1',
                '3de2',
                '3de3',
                '3de4',
                '3de5',
                '3de6',
                '3de7',
                '3p7p',
                '3p7q',
                '3p7r',
                '3p7s',
                '3t6w',
                '3t6x',
                '3t6z',
                '3t71',
                '4cw2',
                '4cw6',
                '4cw3',
                '4d13',
                '4d17',
                '4d19',
                '4h8x',
                '4h8y',
                '4h8z',
                '4h90',
                '4h91',
                '4h92',
                '4h93',
                '4h94',
                '4h9a',
                '4h9b',
                '4h9c',
                '4h9e',
                '4h9f',
                '4h9h',
                '4h9i',
                '4qwm',
                '4ubi',
                '4ubj',
                '4ubk',
                '4ubm',
                '4ubl',
                '4yso',
                '4ysp',
                '4ysq',
                '4ysr',
                '4yss',
                '4yst',
                '4ysu',
                '5gx1',
                '5gx2',
                '5gx3',
                '5gx4',
                '5gx5',
                '5kul',
                '5kun',
                '5kuo',
                '5kuq',
                '5kur',
                '5kus',
                '5kuu',
                '5kuv',
                '5kuw',
                '5kvw',
                '5kvx',
                '5kvz',
                '5kw0',
                '5mcc',
                '5mcd',
                '5mce',
                '5mcf',
                '5mch',
                '5mci',
                '5mcj',
                '5mck',
                '5mcl',
                '5mcm',
                '5mcn',
                '6f7r',
                '6f81',
                '6f82',
                '6nsw',
                '6nsy',
                '6nsz',
                '6nt0',
                '6nt1',
                '6qrr',
                '6qrs',
                '6qrt',
                '6qru',
                '6qrv',
                '6qrw',
                '6qrx',
                '6rgh',
                '6rgp',
                '6rhh',
                '6rhi',
                '6rho',
                '6rhr',
                '6rhu',
                '6rhx',
                '6ri0',
                '6ri2',
                '6ri4',
                '6ri6',
                '6ri8',
                '6rii',
                '6rik']

    return pdb_list


def check_rsrz_outliers(pdb_code):
    """
    Finds the Wilson B value as recorded in the full PDB report
    """

    wilson_b = np.nan

    pdb_info = requests.get(
        'http://files.rcsb.org/pub/pdb/validation_reports/{}/{}/{}_validation.xml.gz'.format(
        pdb_code.lower()[1:3], pdb_code.lower(), pdb_code.lower())
    )
    if pdb_info.status_code == 200:
        try:
            wilson_b = float(pdb_info.text.split('" WilsonBestimate=')[1].split('"')[1])
        except (ValueError, IndexError):
            pass

    return wilson_b


def find_deposition_year(cif_lines):
    """
    Finds year of deposition in mmCIF file header
    """

    year = np.nan

    for subsection in cif_lines:
        if '\n_pdbx_database_status.recvd_initial_deposition_date ' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_pdbx_database_status.recvd_initial_deposition_date '):
                    try:
                        date = line.replace('_pdbx_database_status.recvd_initial_deposition_date ', '')
                        year = int(date.split('-')[0])
                    except ValueError:
                        pass
    return year


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
        if '\n_refine.ls_d_res_high ' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_d_res_high '):
                    try:
                        resolution = float(line.replace('_refine.ls_d_res_high ', ''))
                    except ValueError:
                        pass
                    break

        # Finds structure Rwork
        if '\n_refine.ls_R_factor_R_work ' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_R_factor_R_work '):
                    try:
                        rwork = float(line.replace('_refine.ls_R_factor_R_work ', ''))
                    except ValueError:
                        pass
                    break

        # Finds structure Rfree
        if '\n_refine.ls_R_factor_R_free ' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_refine.ls_R_factor_R_free '):
                    try:
                        rfree = float(line.replace('_refine.ls_R_factor_R_free ', ''))
                    except ValueError:
                        pass
                    break

        # Finds temperature of data collection
        if '\n_diffrn.ambient_temp ' in subsection:
            sub_lines = subsection.split('\n')
            for line in sub_lines:
                if line.startswith('_diffrn.ambient_temp '):
                    try:
                        temperature = float(line.replace('_diffrn.ambient_temp ', ''))
                    except ValueError:
                        pass
                    break

        # Finds residues in asymmetric unit
        if '\n_pdbx_poly_seq_scheme.' in subsection:
            sub_lines = [line for line in subsection.split('\n')
                         if not line.strip() in ['', 'loop_']]
            prop_indices = {}
            prop_num = 0
            for line in sub_lines:
                if line.startswith('_pdbx_poly_seq_scheme.'):
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
        if '\n_struct_conn.' in subsection:
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


def find_pdb_metadata():
    """
    Writes properties of selected protein stuctures to a dataframe ()
    """

    pdb_list = gen_pdb_list()

    pdb_dir = '/home/shared/structural_bioinformatics/rcsb_pdb/mmCIF'
    work_dir = '/home/shared/structural_bioinformatics/RABDAM_PDB-REDO_analysis/PDB_parsing_output/'

    structure_count = 0
    for pdb_code in pdb_list:
        # Tracks progress
        structure_count += 1
        print('{}: {}'.format(pdb_code, (structure_count/len(pdb_list))))

        if os.path.isfile('{}/PDB_rad_dam_series_metadata.pkl'.format(work_dir)):
            all_pdbs_df = pd.read_pickle('{}/PDB_rad_dam_series_metadata.pkl'.format(work_dir))
        else:
            all_pdbs_df = pd.DataFrame({'PDB code': [],
                                        'Deposition year': [],
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
                                        'Per-atom refined B-factors': []})

        # Opens mmCIF PDB file
        pdb_lines = None
        try:
            with open(
                '{}/{}/{}.cif'.format(
                pdb_dir, pdb_code[1:3].lower(), pdb_code.lower()
                ), 'r'
            ) as f:
                pdb_lines = f.read().split('#')
        except FileNotFoundError:
            raise FileNotFoundError('File {}/{}/{}.cif doesn\'t exist'.format(
                pdb_dir, pdb_code[1:3].lower(), pdb_code.lower()
            ))

        # Generates dictionary of properties for each atom
        atom_lines = []
        for subsection in pdb_lines:
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
                            if atom_props['type_symbol'] != 'H':
                                atom_lines.append(atom_props)

        # Finds structural properties of interest
        wilson_b = check_rsrz_outliers(pdb_code)
        year = find_deposition_year(pdb_lines)
        (resolution, rwork, rfree, temperature, seqres, ss_bonds
        ) = find_structurewide_properties(pdb_lines)
        (size, glu_asp_count, glu_asp_ld_count, glu_asp_percent,
         non_canonical_aa_count
        ) = find_aa_properties(seqres, mass_dict)
        (num_terminal_o_atoms, sub_1_any, sub_1_asp_glu, sub_1_disulfide,
         single_model) = find_peratom_properties(atom_lines, seqres, ss_bonds)
        per_atom_b_factors = find_non_per_atom_b_factors(atom_lines, seqres)

        # Summarises PDB features in dataframe
        print('Summarising info for {}'.format(pdb_code))
        indv_pdb_df = pd.DataFrame({'PDB code': [pdb_code],
                                    'Deposition year': [year],
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
                                    'Per-atom refined B-factors': [per_atom_b_factors]})

        all_pdbs_df = pd.concat([all_pdbs_df, indv_pdb_df], axis=0, ignore_index=True)
        all_pdbs_df.to_pickle('{}/PDB_rad_dam_series_metadata.pkl'.format(work_dir))
        all_pdbs_df.to_excel('{}/PDB_rad_dam_series_metadata.xlsx'.format(work_dir), index=False)
