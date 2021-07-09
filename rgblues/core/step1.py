# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Perform EDR3 query
"""

import sys

from astropy import units as u
from astropy.table import Column
from astroquery.gaia import Gaia
import numpy as np


def step1(ra_center, dec_center, search_radius, g_limit, verbose):
    """Perform EDR3 query

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
    r_edr3 : astropy Table
        Table containing the query result.
    nstars : int
        Number of stars in 'r_edr3'.
    nstars_colorcut: int
        Number of stars outside the colour cut
        -0.5 < G_BP - G_RP < 2.0 mag.

    """
    query = f"""
    SELECT source_id, ra, dec,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag

    FROM gaiaedr3.gaia_source
    WHERE 1=CONTAINS(
      POINT('ICRS', {ra_center}, {dec_center}), 
      CIRCLE('ICRS',ra, dec, {search_radius}))
    AND phot_g_mean_mag IS NOT NULL 
    AND phot_bp_mean_mag IS NOT NULL 
    AND phot_rp_mean_mag IS NOT NULL
    AND phot_g_mean_mag < {g_limit}

    ORDER BY ra
    """
    sys.stdout.write('<STEP1> Starting cone search in Gaia EDR3... (please wait)\n  ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_edr3 = job.get_results()
    # compute G_BP - G_RP colour
    r_edr3.add_column(
        Column(r_edr3['phot_bp_mean_mag'] - r_edr3['phot_rp_mean_mag'],
               name='bp_rp', unit=u.mag)
    )
    # colour cut in BP-RP
    mask_colour = np.logical_or((r_edr3['bp_rp'] <= -0.5), (r_edr3['bp_rp'] >= 2.0))
    r_edr3_colorcut = r_edr3[mask_colour]
    nstars = len(r_edr3)
    print(f'        --> {nstars} stars found')
    nstars_colorcut = len(r_edr3_colorcut)
    print(f'        --> {nstars_colorcut} stars outside -0.5 < G_BP-G_RP < 2.0')
    if nstars == 0:
        raise SystemExit('ERROR: no stars found. Change search parameters!')
    if verbose:
        r_edr3.pprint(max_width=1000)

    return r_edr3, nstars, nstars_colorcut
