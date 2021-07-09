# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Compute intersection of EDR3 query with the 15M star sample
"""

import sys

import numpy as np


def step3(r_edr3, edr3_source_id_15M_allsky, debug):
    """Compute intersection of EDR3 query with the 15M star sample

    Parameters
    ----------
    r_edr3 : astropy Table
        Table containing the initial EDR3 query result.
    edr3_source_id_15M_allsky : numpy array
        Array with source_id values for the 15M star sample.
    debug : bool
        If True, display debugging information.

    Returns
    -------
    intersection : numpy array
        Array with the source_id values of the intersection between
        the EDR3 query and the 15M star sample.

    """
    sys.stdout.write('<STEP3> Cross-matching EDR3 with 15M subsample... (please wait)')
    sys.stdout.flush()
    set1 = set(np.array(r_edr3['source_id']))
    set2 = set(edr3_source_id_15M_allsky)
    intersection = set2.intersection(set1)
    print(f'\n        --> {len(intersection)} stars in common with 15M sample')
    if debug:
        print(len(set1), len(set2), len(intersection))

    return intersection
