# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Retrieve additional parameters from the StarHorse sample
"""

from astropy.table import vstack, join
import pyvo
import requests


def step2(r_edr3, starhorse_block, verbose, debug):
    """Retrieve additional parameters from the StarHorse sample

    Parameters
    ----------
    r_edr3 : astropy Table
        Table containing the initial EDR3 query result.
    starhorse_block : int
        Integer number indicating the maximum number of stars whose
        parameters are retrieved in each single query to Gaia@AIP.
        This number should not be too large to avoid an out-of-time
        error (a typical useful value is 100).
    verbose : bool
        If True, display additional information.
    debug : bool
        If True, display the intermediate tables retrieved in the
        different queries to Gaia@AIP.

    Returns
    -------
    r_edr3 : astropy Table
        Updated table containing the additional StarHorse parameters.

    """
    param_starhorse = [
        'dr3_source_id', 'sh_gaiaflag', 'sh_outflag',
        'dist05', 'dist16', 'dist50', 'dist84', 'dist95',
        'av05', 'av16', 'av50', 'av84', 'av95',
        'teff16', 'teff50', 'teff84',
        'logg16', 'logg50', 'logg84',
        'met16', 'met50', 'met84',
        'mass16', 'mass50', 'mass84',
        'xgal', 'ygal', 'zgal', 'rgal',
        'ruwe', 'angular_distance', 'magnitude_difference',
        'proper_motion_propagation', 'dup_max_number'
    ]
    print('<STEP2> Retrieving StarHorse data from Gaia@AIP... (please wait)')
    print(f'        pyvo version {pyvo.__version__}')
    print(f'        TAP service GAIA@AIP')
    nstars_per_block = starhorse_block
    nstars = len(r_edr3)
    nblocks = int(nstars / nstars_per_block)
    r_starhorse = None
    if nstars - nblocks * nstars_per_block > 0:
        nblocks += 1
    for iblock in range(nblocks):
        irow1 = iblock * nstars_per_block
        irow2 = min(irow1 + nstars_per_block, nstars)
        print(f'        Starting query #{iblock + 1} of {nblocks}...')
        dumstr = ','.join([str(item) for item in r_edr3[irow1:irow2]['source_id']])
        query = f"""
        SELECT {','.join(param_starhorse)}
        FROM gaiadr2_contrib.starhorse
        WHERE dr3_source_id IN ({dumstr})
        """
        tap_session = requests.Session()
        tap_session.headers['Authorization'] = "kkk"
        tap_service = pyvo.dal.TAPService('https://gaia.aip.de/tap', session=tap_session)
        tap_result = tap_service.run_sync(query)
        if debug:
            print(tap_result.to_table())
        if iblock == 0:
            r_starhorse = tap_result.to_table()
        else:
            r_starhorse = vstack([r_starhorse, tap_result.to_table()],
                                 join_type='exact', metadata_conflicts='silent')

    nstars_starhorse = len(r_starhorse)
    if verbose:
        if nstars_starhorse > 0:
            r_starhorse.pprint(max_width=1000)

    # join tables
    print(f'        --> {nstars_starhorse} stars found in StarHorse')
    print('        Joining EDR3 and StarHorse queries...')
    r_starhorse.rename_column('dr3_source_id', 'source_id')
    r_edr3 = join(r_edr3, r_starhorse, keys='source_id', join_type='outer')
    r_edr3.sort('ra')
    if verbose:
        r_edr3.pprint(max_width=1000)

    return r_edr3
