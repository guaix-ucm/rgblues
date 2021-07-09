# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Generate output CSV files
"""

import glob
import os
import sys

from astropy.table import Column
import numpy as np

OUTTYPES = ['edr3', '15m', 'var']


def step7(r_edr3, basename, starhorse_block, nvariables, r_cross_var,
          intersection, table15M, verbose, debug):
    """Generate output CSV files

    Parameters
    ----------
    r_edr3 : astropy Table
        Table containing the EDR3 query.
    basename : str
        Base name for output files.
    starhorse_block : int
        Integer number indicating the maximum number of stars whose
        parameters are retrieved in each single query to Gaia@AIP.
        This number should not be too large to avoid an out-of-time
        error (a typical useful value is 100).
    nvariables : int
        Final number of variable stars in EDR3 query.
    r_cross_var : astropy Table
        Table with the cross-match of the DR2 and EDR3 queries.
    intersection : numpy array
        Array with the source_id values of the intersection between
        the EDR3 query and the 15M star sample.
    table15M : astropy FITS binary table
        Table containing parameters for the stars in the 15M sample.
    verbose : bool
        If True, display additional information.
    debug : bool
        If True, display debugging information.

    """
    sys.stdout.write('<STEP7> Saving output CSV files...')
    sys.stdout.flush()
    r_edr3.add_column(Column(np.zeros(len(r_edr3)), name='number_csv', dtype=int))
    for item in OUTTYPES:
        r_edr3.add_column(Column(np.zeros(len(r_edr3)), name=f'number_{item}', dtype=int))
    outlist = [f'./{basename}_{ftype}.csv' for ftype in OUTTYPES]
    filelist = glob.glob('./*.csv')
    # remove previous versions of the output files (if present)
    for file in outlist:
        if file in filelist:
            try:
                os.remove(file)
            except OSError:
                print(f'ERROR: while deleting existing file {file}')
    # columns to be saved (use a list to guarantee the same order)
    outcolumns_list = ['source_id', 'ra', 'dec', 'b_rgb', 'g_rgb', 'r_rgb', 'g_br_rgb',
                       'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag']
    # define column format with a dictionary
    outcolumns = {'source_id': '19d',
                  'ra': '14.9f',
                  'dec': '14.9f',
                  'b_rgb': '6.2f',
                  'g_rgb': '6.2f',
                  'r_rgb': '6.2f',
                  'g_br_rgb': '6.2f',
                  'phot_g_mean_mag': '8.4f',
                  'phot_bp_mean_mag': '8.4f',
                  'phot_rp_mean_mag': '8.4f'
                  }
    if set(outcolumns_list) != set(outcolumns.keys()):
        raise SystemExit('ERROR: check outcolumns_list and outcolumns')
    csv_header_ini = 'number,' + ','.join(outcolumns_list)
    flist = []
    for ftype in OUTTYPES:
        f = open(f'{basename}_{ftype}.csv', 'wt')
        flist.append(f)
        if (starhorse_block > 0) and (ftype in ['edr3', '15m']):
            if debug and (ftype == '15m'):
                csv_header = csv_header_ini + \
                             ',av50,met50,dist50,b_rgb_bis,g_rgb_bis,r_rgb_bis,g_br_rgb_bis,' \
                             'phot_g_mean_mag_bis,phot_bp_mean_mag_bis,phot_rp_mean_mag_bis,' \
                             'av50_bis,met50_bis,dist50_bis'
            else:
                csv_header = csv_header_ini + ',av50,met50,dist50'
        else:
            csv_header = csv_header_ini
        f.write(csv_header + '\n')
    # save each star in its corresponding output file
    krow = np.ones(len(OUTTYPES), dtype=int)
    for irow, row in enumerate(r_edr3):
        cout = []
        for item in outcolumns_list:
            cout.append(eval("f'{row[item]:" + f'{outcolumns[item]}' + "}'"))
        iout = 0
        if nvariables > 0:
            if row['source_id'] in r_cross_var['dr3_source_id']:
                iout = 2
        if iout == 0:
            if starhorse_block > 0:
                for item in ['av50', 'met50', 'dist50']:
                    value = row[item]
                    if isinstance(value, float):
                        pass
                    else:
                        value = 99.999
                    cout.append(f'{value:7.3f}')
            if row['source_id'] in intersection:
                iout = 1
                if debug:
                    iloc = np.argwhere(table15M['source_id'] == row['source_id'])[0][0]
                    cout.append(f"{table15M['B_rgb'][iloc]:6.2f}")
                    cout.append(f"{table15M['G_rgb'][iloc]:6.2f}")
                    cout.append(f"{table15M['R_rgb'][iloc]:6.2f}")
                    cout.append(f"{table15M['G_BR_rgb'][iloc]:6.2f}")
                    cout.append(f"{table15M['G_gaia'][iloc]:8.4f}")
                    cout.append(f"{table15M['BP_gaia'][iloc]:8.4f}")
                    cout.append(f"{table15M['RP_gaia'][iloc]:8.4f}")
                    cout.append(f"{table15M['av50'][iloc]:7.3f}")
                    cout.append(f"{table15M['met50'][iloc]:7.3f}")
                    cout.append(f"{table15M['dist50'][iloc]:7.3f}")
        flist[iout].write(f'{krow[iout]:6d}, ' + ','.join(cout) + '\n')
        r_edr3[irow]['number_csv'] = iout
        r_edr3[irow][f'number_{OUTTYPES[iout]}'] = krow[iout]
        krow[iout] += 1
    for f in flist:
        f.close()
    print('OK')

    if verbose:
        print(r_edr3)
