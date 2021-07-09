# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
DR2 query to identify variable stars
"""

import sys

from astroquery.gaia import Gaia


def step4(ra_center, dec_center, search_radius, g_limit, verbose):
    """Perform DR2 query

    Parameters
    ----------
    ra_center : float
        Right ascension (decimal degree) corresponding to the center
        of the field of view.
    dec_center : float
        Declination (decimal degree) corresponding to the center
        of the field of view.
    search_radius : float
        Radius (decimal degrees) of the field of view.
    g_limit : float
        Limiting Gaia G magnitude.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_dr2 : astropy Table
        Table containing the query result.
    nvariables : int
        Number of variable stars found.
    mask_var : numpy array
        Array indicating whether the star is variable or not.

    """
    query = f"""
    SELECT source_id, ra, dec, phot_g_mean_mag, phot_variable_flag

    FROM gaiadr2.gaia_source
    WHERE  1=CONTAINS(
      POINT('ICRS', {ra_center}, {dec_center}), 
      CIRCLE('ICRS',ra, dec, {search_radius}))
    AND phot_g_mean_mag < {g_limit}
    """
    sys.stdout.write('<STEP4> Looking for variable stars in Gaia DR2... (please wait)\n  ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_dr2 = job.get_results()
    nstars_dr2 = len(r_dr2)
    if nstars_dr2 == 0:
        nvariables = 0
        mask_var = None
    else:
        if isinstance(r_dr2['phot_variable_flag'][0], bytes):
            mask_var = r_dr2['phot_variable_flag'] == b'VARIABLE'
        elif isinstance(r_dr2['phot_variable_flag'][0], str):
            mask_var = r_dr2['phot_variable_flag'] == 'VARIABLE'
        else:
            raise SystemExit('Unexpected type of data in column phot_variable_flag')
        nvariables = sum(mask_var)
        print(f'        --> {nstars_dr2} stars in DR2, ({nvariables} initial variables)')
    if nvariables > 0:
        if verbose:
            r_dr2[mask_var].pprint(max_width=1000)

    return r_dr2, nvariables, mask_var
