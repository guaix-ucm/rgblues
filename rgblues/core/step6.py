# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Compute RGB magnitudes
"""

import sys

from astropy import units as u
from astropy.table import Column
import numpy as np
from numpy.polynomial import Polynomial


def step6(r_edr3, verbose):
    """Cross-matching between DR2 and EDR3 queries

    Parameters
    ----------
    r_edr3 : astropy Table
        Table containing the EDR3 query.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_edr3 : astropy Table
        Updated EDR3 table including the predicted RGB magnitudes.

    """
    sys.stdout.write('<STEP6> Computing RGB magnitudes...')
    sys.stdout.flush()
    # predict RGB magnitudes
    coef_B = np.array([-0.13748689, 0.44265552, 0.37878846, -0.14923841, 0.09172474, -0.02594726])
    coef_G = np.array([-0.02330159, 0.12884074, 0.22149167, -0.1455048, 0.10635149, -0.0236399])
    coef_R = np.array([0.10979647, -0.14579334, 0.10747392, -0.1063592, 0.08494556, -0.01368962])
    coef_X = np.array([-0.01252185, 0.13983574, 0.23688188, -0.10175532, 0.07401939, -0.0182115])

    poly_B = Polynomial(coef_B)
    poly_G = Polynomial(coef_G)
    poly_R = Polynomial(coef_R)
    poly_X = Polynomial(coef_X)

    r_edr3.add_column(
        Column(np.round(r_edr3['phot_g_mean_mag'] + poly_B(r_edr3['bp_rp']), 2),
               name='b_rgb', unit=u.mag, format='.2f'), index=3
    )
    r_edr3.add_column(
        Column(np.round(r_edr3['phot_g_mean_mag'] + poly_G(r_edr3['bp_rp']), 2),
               name='g_rgb', unit=u.mag, format='.2f'), index=4
    )
    r_edr3.add_column(
        Column(np.round(r_edr3['phot_g_mean_mag'] + poly_R(r_edr3['bp_rp']), 2),
               name='r_rgb', unit=u.mag, format='.2f'), index=5
    )
    r_edr3.add_column(
        Column(np.round(r_edr3['phot_g_mean_mag'] + poly_X(r_edr3['bp_rp']), 2),
               name='g_br_rgb', unit=u.mag, format='.2f'), index=6
    )
    print('OK')
    if verbose:
        r_edr3.pprint(max_width=1000)

    return r_edr3
