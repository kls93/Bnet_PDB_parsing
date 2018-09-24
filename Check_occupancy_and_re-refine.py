
import os
import requests
import shutil
import numpy as np
import pandas as pd

refmac_input = ('make check NONE\n'
                'make -\n'
                '    hydrogen ALL -\n'
                '    hout NO -\n'
                '    peptide NO -\n'
                '    cispeptide YES -\n'
                '    ssbridge YES -\n'
                '    symmetry YES -\n'
                '    sugar YES -\n'
                '    connectivity NO -\n'
                '    link NO\n'
                'refi -\n'
                '    type REST -\n'
                '    resi MLKF -\n'
                '    meth CGMAT -\n'
                '    bref ISOT\n'
                'ncyc 10\n'
                'scal -\n'
                '    type SIMP -\n'
                '    LSSC -\n'
                '    ANISO -\n'
                '    EXPE\n'
                'solvent YES\n'
                'weight -\n'
                '    AUTO\n'
                'monitor MEDIUM -\n'
                '    torsion 10.0 -\n'
                '    distance 10.0 -\n'
                '    angle 10.0 -\n'
                '    plane 10.0 -\n'
                '    chiral 10.0 -\n'
                '    bfactor 10.0 -\n'
                '    bsphere 10.0 -\n'
                '    rbond 10.0 -\n'
                '    ncsr 10.0\n'
                'labin  FP=FP SIGFP=SIGFP -\n'
                '   FREE=FREE\n'
                'labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT '
                'DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\n'
                'PNAME unknown\n'
                'DNAME unknown\n'
                'RSIZE 80\n'
                'EXTERNAL WEIGHT SCALE 10.0\n'
                'EXTERNAL USE MAIN\n'
                'EXTERNAL DMAX 4.2\n'
                'END')
with open('/Users/ks17361/Downloads/REFMAC_input.txt', 'w') as f:
    f.write(refmac_input)
refmac_input_file = '/Users/ks17361/Downloads/REFMAC_input.txt'
refmac_output_file = '/Users/ks17361/Downloads/REFMAC_output.txt'

rad_dam_pdbs = ['1N6X', '1N6Y', '2BHX', '2BI1', '2BI2', '2BI3', '2BI5', '2BI9',
                '2BIA', '3GHR', '3GHS', '3GHT', '3GHU', '3MNB', '3MNC', '3MNS',
                '3MNX', '3MO3', '3MO6', '3MO9', '3MOC', '3MTY', '3ODF', '3MU0',
                '3MU1', '3P7P', '3P7Q', '3P7R', '3P7S', '3P7T', '3P7U', '3P7V',
                '3P7W', '4EK0', '4EKA', '4EKB', '4EKH', '4EKO', '4EKT', '4EL2',
                '4EL3', '4EL7', '4ELA', '4EP8', '4EPB', '4EPD', '4EPE', '4M4M',
                '4M4F', '4M4H', '4M4I', '4M4J', '4M4L', '4X4B', '4X4C', '4X4D',
                '4X4E', '4X4F', '4X4G', '4X4H', '4X4I', '5EEU', '5EEV', '5EEW',
                '5EEX', '5EEY', '5EEZ', '5EF0', '5EF1', '5EF2', '5EF3', '2J5K',
                '2J5Q', '2J5R', '2BLO', '2BLQ', '2BLP', '2BLZ', '2BLR', '2BLU',
                '2BLV', '2BLW', '2BLX', '2BLY', '2BN3', '2BN1', '1QID', '1QIE',
                '1QIF', '1QIG', '1QIH', '1QII', '1QIJ', '1QIK', '1QIM', '2BXY',
                '2BXZ', '2BY0', '2BY1', '2BY2', '2BY3', '2BY5', '2BY6', '2BY7',
                '2BY8', '2BY9', '2BYA', '3TMW', '3TMX', '3TMU', '3TMV', '4V2L',
                '4V2M', '4V2N', '5MCC', '5MCD', '5MCE', '5MCF', '5MCH', '5MCI',
                '5MCJ', '5MCK', '5MCL', '5MCM', '5MCN', '3ONC', '3ONB', '5KUL',
                '5KUN', '5KUO', '5KUQ', '5KUR', '5KUS', '5KUU', '5KUV', '5KUW',
                '5KUZ', '5KV0', '5KV1', '5KV2', '5KV3', '5KV4', '5KV5', '5KV6',
                '5KV7', '5KVW', '5KVX', '5KVZ', '5KW0', '5KW3', '5KW4', '5KW5',
                '5KW7', '5KW8', '5KXK', '5KXL', '5KXM', '5KXN', '5KXO', '5KXP',
                '5KXR', '5KXS', '5KXT', '5KXW', '5KXX', '5KXY', '5KXZ', '5KY1',
                '3NSB', '3V7V', '3VCH', '3V88', '3V87', '3VCI', '3V8A', '3V82',
                '3VCJ', '3VCE', '3V84', '3VCK', '3VCG', '5L9J', '5LA5', '5LA8',
                '5LAF', '5LAG', '5LAN', '2XBR', '2XBS', '2W1L', '2W1M', '2W1X',
                '2W1Y', '4GCC', '4GCD', '4GCE', '4GCF', '4ETA', '4ETB', '4ETC',
                '4ETD', '4ETE', '4ET8', '4ET9', '5LH0', '5LH1', '5LN0', '5LH3',
                '5LH5', '5LMH', '5LH6', '5LH7', '5M3S', '5O41', '2YBH', '2YBI',
                '2YBJ', '2YBL', '2YBM', '2YBN', '1OWM', '1OWN', '1OWO', '1OWP',
                '2VO2', '2VNX', '2WJM', '2WJN', '2YAE', '2YAF', '2YAH', '2YAM',
                '2YAO', '2YAP', '2YAQ', '2YAR', '3T6W', '3T6X', '3T6Z', '3T71',
                '3VKP', '3VKQ', '3VKR', '4QWM', '4UBI', '4UBJ', '4UBK', '4UBL',
                '4UBM', '4UBN', '4UBO', '4W1P', '4W1Q', '4W1R', '4W1S', '4CW2',
                '4CW3', '4CW6', '4D13', '4D17', '4D19', '4UPV', '4UQL', '1H5D',
                '1H5E', '1H5F', '1H5G', '1H5H', '1H5I', '1H5J', '1H5K', '1H5L',
                '5GX1', '5GX2', '5GX3', '5GX4', '5GX5', '5I6K', '5I6L', '5I6M',
                '5I6N', '5I6O', '5I6P', '4YSO', '4YSP', '4YSQ', '4YSR', '4YSS',
                '4YST', '4YSU', '4H8X', '4H8Y', '4H8Z', '4H90', '4H91', '4H92',
                '4H93', '4H94', '4H9A', '4H9B', '4H9C', '4H9E', '4H9F', '4H9H',
                '4H9I']
xfel_pdbs = ['5wrc', '5b1e', '4n5r', '5b35', '5cnf', '5wr8', '5l8m', '5hdd',
             '4tnl', '4rwd', '4ixq', '4yop', '5trx', '4ym8', '6b6e', '4tnk',
             '5d9c', '5ws6', '5gti', '5nj4', '3pcq', '5oer', '4wl9', '5wrb',
             '5wr9', '5fgt', '5hdc', '4ysc', '5k2a', '3wxu', '5nm4', '5o64',
             '5d9b', '5xfc', '5cn8', '5h2k', '5b1g', '6g7h', '5y5j', '5x1f',
             '3wxq', '6ch7', '5kxu', '5c6j', '5ndc', '5h2l', '5w0p', '4rvy',
             '5fvf', '4w4q', '5oqe', '4wla', '5foy', '5hd5', '5xfe', '6b6f',
             '5br8', '4o9r', '4tnh', '5h2n', '5k2c', '5b1d', '5m7k', '4ub6',
             '5x1b', '5cnc', '5k2b', '6b5y', '5jd2', '5tis', '5ws5', '6b69',
             '5d4i', '4pbu', '5kaf', '3wg7', '5b6x', '4ub8', '5onx', '5mqw',
             '5g0z', '4rw2', '5o8c', '5j7a', '4qx3', '5y5h', '5ooz', '5um1',
             '5oq9', '4qx2', '6b68', '5xjv', '5xfd', '5mg1', '5hl4', '4yay',
             '5cn5', '5cn9', '5mg0', '4ysa', '4q70', '5b6y', '5o8a', '6g7l',
             '5h2i', '5ejx', '4zix', '5y5l', '5o8b', '5wr4', '5unf', '4yup',
             '4et9', '6g7i', '6b6a', '5xez', '5e79', '4hwy', '5y5k', '5gth',
             '4tni', '4fby', '5mnd', '5ws3', '5x2i', '5hqw', '5swe', '4nc3',
             '5f81', '5e54', '5cn7', '5dlh', '4tnj', '4ixr', '6b6d', '6g7k',
             '6apk', '5cn6', '5h2p', '4cas', '6b6b', '5joo', '5wr2', '5m7j',
             '5y5i', '5e7c', '5kj7', '6b5x', '4s1l', '5b34', '5osn', '4qx0',
             '5ung', '4qx1', '5x19', '5jom', '5wr3', '5cmv', '5foz', '5tud',
             '4xou', '5kai', '6b6c', '4z8k', '5ttc', '4zwj', '4ziz', '3wxs',
             '5k2d', '5hqd', '5b6v', '5h2m', '5swd', '3wul', '5b6w', '5u5q',
             '4ac5', '4uyo', '5h2o', '3wum', '5hd3', '5lbr', '5b6z', '5wra',
             '5g37', '5cnd', '5cne', '5dm9', '5v56', '5gyz', '5o89', '5c6l',
             '3wun', '4rw1', '4et8', '5cn4', '5ccg', '5oqa', '4zqx', '6g7j',
             '5h2h', '5w97', '5hds', '4ow3', '5c6i', '3wxt', '5cnb', '5b1f',
             '5y5m', '4pnj', '5h2j', '5o4c', '5cng']

pdb_lists = {'Rad_dam_series': rad_dam_pdbs,
             'XFEL': xfel_pdbs}

for dir_name, pdb_list in pdb_lists.items():
    if os.path.isdir('/Users/ks17361/Downloads/{}'.format(dir_name)):
        shutil.rmtree('/Users/ks17361/Downloads/{}'.format(dir_name))
    os.mkdir('/Users/ks17361/Downloads/{}'.format(dir_name))
    os.mkdir('/Users/ks17361/Downloads/{}/Reprocessed'.format(dir_name))

    raw_pdbs = []
    reprocessed_pdbs = []
    unprocessable_pdbs = []

    for pdb in pdb_list:
        pdb = pdb.lower()
        print(pdb)

        pdb_file_lines = requests.get(
            'http://www.rcsb.org/pdb/files/{}.pdb'.format(pdb)
        ).text.split('\n')

        ss_bonds = []
        new_pdb_lines = []
        reprocessed = False

        for line in pdb_file_lines:
            if line[0:6] == 'SSBOND':
                ss_bonds.append(line[15:22].replace(' ', ''))
                ss_bonds.append(line[29:36].replace(' ', ''))

            if (
                    line[0:4] == 'ATOM'
                and line[17:20] in ['ASP', 'GLU']
                and float(line[54:60]) != 1.0
                and line[16:17] == ' '
            ):
                reprocessed = True

                new_pdb_line = line[:54] + '  1.00' + line[60:]
                new_pdb_lines.append(new_pdb_line)
            elif (
                    line[0:4] == 'ATOM'
                and line[17:20] == 'CYS'
                and line[21:27].replace(' ', '') in ss_bonds
                and float(line[54:60]) != 1.0
                and line[16:17] == ' '
            ):
                reprocessed = True

                new_pdb_line = line[:54] + '  1.00' + line[60:]
                new_pdb_lines.append(new_pdb_line)
            else:
                new_pdb_lines.append(line)

        atom_ids = []
        conformers = [[] for i in range(len(new_pdb_lines))]
        occupancy = [[] for i in range(len(new_pdb_lines))]

        for line in new_pdb_lines:
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
                if not atom_id in atom_ids:
                    atom_ids.append(atom_id)
                index = atom_ids.index(atom_id)
                if not line[16:17] in conformers[index]:
                    conformers[index].append(line[16:17])
                    occupancy[index].append(float(line[54:60]))

        unprocessable = False
        for index, atom_id in enumerate(atom_ids):
            if np.sum(occupancy[index]) != 1.0:
                unprocessable = True
            elif (
                np.sum(occupancy[index]) == 1.0
                and atom_id.split('_')[0] in ss_bonds
            ):
                unprocessable = True

        if reprocessed is False and unprocessable is False:
            raw_pdbs.append(pdb)
        elif reprocessed is True and unprocessable is False:
            reprocessed_pdbs.append(pdb)
            with open('/Users/ks17361/Downloads/{}/Reprocessed/{}.pdb'.format(dir_name, pdb), 'w') as f:
                for line in new_pdb_lines:
                    f.write('{}\n'.format(line))
        else:
            unprocessable_pdbs.append(pdb)

    for pdb in reprocessed_pdbs:
        xyzin = '/Users/ks17361/Downloads/{}/Reprocessed/{}.pdb'.format(dir_name, pdb)
        xyzout = '/Users/ks17361/Downloads/{}/Reprocessed/{}_refmac.pdb'.format(dir_name, pdb)

        hklin = '/Users/ks17361/Downloads/{}/Reprocessed/{}.mtz'.format(dir_name, pdb)
        os.system('curl http://www.cmbi.ru.nl/pdb_redo/{}/{}/{}_final.mtz --output {}'.format(
           pdb[1:3], pdb, pdb, hklin
        ))
        hklout = '/Users/ks17361/Downloads/{}/Reprocessed/{}_refmac.mtz'.format(dir_name, pdb)

        if requests.get('http://www.cmbi.ru.nl/pdb_redo/{}/{}/{}_final.mtz'.format(pdb[1:3], pdb, pdb)).status_code < 300:
            print('Re-refining {}'.format(pdb))
            os.system('refmac5 XYZIN {} XYZOUT {} HKLIN {} HKLOUT {} < {} > {}'.format(
                xyzin, xyzout, hklin, hklout, refmac_input_file, refmac_output_file
            ))
        else:
            unprocessable_pdbs.append(pdb)

    with open('/Users/ks17361/Downloads/{}/PDBs_no_reprocessing.txt'.format(dir_name), 'w') as f:
        for pdb in raw_pdbs:
            f.write('{}\n'.format(pdb))
    with open('/Users/ks17361/Downloads/{}/Automatically_reprocessed_PDBs.txt'.format(dir_name), 'w') as f:
        for pdb in reprocessed_pdbs:
            f.write('{}\n'.format(pdb))
    with open('/Users/ks17361/Downloads/{}/PDBs_for_manual_reprocessing.txt'.format(dir_name), 'w') as f:
        for pdb in unprocessable_pdbs:
            f.write('{}\n'.format(pdb))

    with open('/Users/ks17361/Downloads/{}/Reprocessed/RABDAM_input.txt'.format(dir_name), 'w') as f:
        for pdb in reprocessed_pdbs:
            f.write('/Users/ks17361/Downloads/{}/Reprocessed/{}_refmac.pdb,\n'.format(dir_name, pdb))

os.remove(refmac_input_file)
os.remove(refmac_output_file)
