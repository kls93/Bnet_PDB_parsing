
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
        except ValueError:
            pass
        try:
            rsrz_outliers_sim_res = float(
                pdb_info.split('" relative-percentile-percent-RSRZ-outliers=')[1].split('"')[1]
            )
        except ValueError:
            pass
        try:
            wilson_b = float(pdb_info.split('" WilsonBestimate=')[1].split('"')[1])
        except ValueError:
            pass

    return rsrz_outliers_all, rsrz_outliers_sim_res, wilson_b


def calc_bnet(pdb_code):
    """
    Runs RABDAM to calculate Bnet value
    """

    bnet = np.nan

    os.system('python /Volumes/Seagate_Backup_Plus_Drive/RABDAM_PDB_run/RABDAM/rabdam/rabdam.py -f {} -o bnet'.format(pdb_code))

    return bnet


def find_structurewide_properties(remark_lines):
    """
    Finds structure properties in the PDB header
    """

    resolution = np.nan
    rwork = np.nan
    rfree = np.nan
    temperature = np.nan
    seqres = {}
    single_model = True

    for line in remark_lines:
        no_whtspc_line = line.replace(' ', '')

        # Checks there is only a single model in the PDB file
        if line[0:6] == 'MODEL ':
            single_model = False

        # Finds structure resolution
        elif (
                no_whtspc_line.startswith('REMARK2')
            and 'RESOLUTION' in no_whtspc_line
            and 'ANGSTROM' in no_whtspc_line
        ):
            try:
                resolution = float(line[23:30])
            except ValueError:
                pass

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
                pass

        # Finds structure Rfree
        elif (
            any (x in no_whtspc_line for x in
                 ['REMARK3FREERVALUE:', 'REMARK3FREERVALUE(NOCUTOFF):'])
        ):
            sub_lines = no_whtspc_line.split(':')
            try:
                rfree = float(sub_lines[-1])
            except ValueError:
                pass

        # Finds temperature of data collection
        elif no_whtspc_line.startswith('REMARK200TEMPERATURE(KELVIN):'):
            sub_lines = no_whtspc_line.split(':')
            try:
                temperature = float(sub_lines[-1])
            except ValueError:
                pass

        # Finds residues in asymmetric unit
        elif line[0:6] == 'SEQRES':
            if not line[11:12] in seqres:
                seqres[line[11:12]] = []
            seqres[line[11:12]] += [res.upper() for res in line[19:].split()]

    return (resolution, rwork, rfree, temperature, seqres, single_model)


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

    num_terminal_o_atoms = 0
    all_reprocessed = False
    asp_glu_reprocessed = False

    alternate_conformers = {}

    for line in atom_lines:
        # Checks ATOM records for sub-1 occupancy non-alternate conformers
        if (
                line[17:20] in [res for chain in seqres.values() for res in chain]
            and float(line[54:60]) != 1.0
        ):
            if line[16:17] == ' ':
                all_reprocessed = True
                if line[17:20] in ['ASP', 'GLU']:
                    asp_glu_reprocessed = True
            elif line[16:17] != ' ':
                atom_id = '_'.join([line[21:27].replace(' ', ''), line[17:20],
                                    line[12:16].strip()])
                if not atom_id in alternate_conformers:
                    alternate_conformers[atom_id] = []
                alternate_conformers[atom_id].append(float(line[54:60]))

        # Finds number of terminal O atoms
        if (
                line[17:20] in ['ASP', 'GLU']
            and line[12:16].strip() in ['OD1', 'OD2', 'OE1', 'OE2']
        ):
            num_terminal_o_atoms += 1

    # Checks that the occupancies of all alternate conformers sum to 1
    for atom_id in alternate_conformers.keys():
        if sum(alternate_conformers[atom_id]) != 1.0:
            all_reprocessed = True
            if atom_id.split('_')[1] in ['ASP', 'GLU']:
                asp_glu_reprocessed = True

    return num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed


def find_non_per_atom_b_factors(atom_lines, seqres):
    """
    Identifies structures with 20% or more residues with N, CA, C and O atoms
    with the same B-factor values
    """

    res_bfactors = {}
    for line in atom_lines:
        if line[17:20] in [res for chain in seqres.values() for res in chain]:
            res_id = line[21:27].replace(' ', '') + '_' + line[16:20]
            atom_id = line[12:16].strip()
            if not res_id in res_bfactors:
                res_bfactors[res_id] = {}
            res_bfactors[res_id][atom_id] = float(line[60:66])

    per_res_bfactor_count = 0
    for res_id in res_bfactors:
        main_chain_bfactors = []

        res_name = res_id[-3:]
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

    with open('../Protein_and_NA_PDB_IDs/Protein_PDB_IDs.txt', 'r') as f:
        protein_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    with open('../Protein_and_NA_PDB_IDs/NA_PDB_IDs.txt', 'r') as f:
        na_pdbs = [pdb.strip('\n') for pdb in f.readlines() if len(pdb) == 5]

    with open('../Protein_and_NA_PDB_IDs/XFEL_PDB_IDs.txt', 'r') as f:
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

    with open('../Progress_track.txt', 'w') as f:
        f.write('PDBs processed so far:\n')
    with open('../Unprocessed_PDBs.txt', 'w') as f:
        f.write('PDBs unable to be processed:\n')

    for pdb_code in protein_pdbs:
        # Marks progress
        structure_count += 1
        print('{}: {}'.format(pdb_code, (structure_count/len(protein_pdbs))))
        with open('../Progress_track.txt', 'w') as f:
            f.write('{}\n'.format(pdb_code))

        # Downloads PDB file of asymmetric unit and extracts
        pdb_request = requests.get(
            'https://files.rcsb.org/download/{}.pdb'.format(pdb_code.lower())
        )
        if pdb_request.status_code != 200:
            with open('../Unprocessed_PDBs.txt', 'w') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to process {}'.format(pdb_code))
            continue

        pdb_lines = pdb_request.text.split('\n')
        atom_lines = []
        remark_lines = []
        remark_end = False
        for line in pdb_lines:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                remark_end = True
                atom_lines.append(line)
            if remark_end is False:
                remark_lines.append(line)

        # Finds structure properties of interest
        electron_density = check_e_dens(pdb_code)
        (rsrz_outliers_all, rsrz_outliers_sim_res, wilson_b
        ) = check_rsrz_outliers(pdb_code)
        bnet = calc_bnet(pdb_code)
        (resolution, rwork, rfree, temperature, seqres, single_model
        ) = find_structurewide_properties(remark_lines)
        if single_model is False:
            with open('../Unprocessed_PDBs.txt', 'w') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: More than one model present in {}'.format(pdb_code))
            continue
        (size, glu_asp_count, glu_asp_percent, non_canonical_aa_count
        ) = find_aa_properties(seqres, mass_dict)
        (num_terminal_o_atoms, all_reprocessed, asp_glu_reprocessed
        ) = find_peratom_properties(atom_lines, seqres)
        per_atom_b_factors = find_non_per_atom_b_factors(atom_lines, seqres)

        # Summarises PDB features in csv file
        if (
               np.isnan(resolution)
            or np.isnan(rwork)
            or np.isnan(rfree)
            or np.isnan(temperature)
            or size == 0
            or num_terminal_o_atoms == 0
            or np.isnan(bnet)
        ):
            with open('../Unprocessed_PDBs.txt', 'w') as f:
                f.write('{}\n'.format(pdb_code))
            print('WARNING: Failed to process {}'.format(pdb_code))
        else:
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
            all_pdbs_df.to_pickle('PDB_file_properties.pkl')

    return all_pdbs_df


def calc_bnet_percentile(all_pdbs_df):
    """
    Calculates bnet_percentile for PDB structure XXXX by comparing its Bnet
    value to the Bnet values of structures of a similar resolution
    """

    for pdb_code in all_pdbs_df['PDB code'].tolist():
        bnet_percentile = np.nan

    all_pdbs_df.to_pickle('PDB_file_properties.pkl')

    return all_pdbs_df
