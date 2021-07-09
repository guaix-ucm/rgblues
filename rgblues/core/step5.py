# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Cross-matching between DR2 and EDR3 to identify the variable stars
"""

import sys

from astroquery.gaia import Gaia


def step5(r_dr2, mask_var, r_edr3, verbose):
    """Cross-matching between DR2 and EDR3 queries

    Parameters
    ----------
    r_dr2 : astropy Table
        Table containing the DR2 query.
    mask_var : numpy array
        Array indicating whether the star is variable or not
        in the DR2 query.
    r_edr3 : astropy Table
        Table containing the EDR3 query.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_cross_var : astropy Table
        Table with the cross-match of the DR2 and EDR3 queries.
    nvariables : int
        Final number of variable stars in EDR3 query.

    """
    # generate sequence of source_id of variable stars
    dumstr = ','.join([str(item) for item in r_dr2[mask_var]['source_id']])
    # cross-match
    query = f"""
    SELECT *
    FROM gaiaedr3.dr2_neighbourhood
    WHERE dr2_source_id IN ({dumstr})
    ORDER BY angular_distance
    """
    sys.stdout.write('<STEP5> Cross-matching variables in DR2 with stars in EDR3... (please wait)\n  ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_cross_var = job.get_results()
    if verbose:
        r_cross_var.pprint(max_width=1000)
    nvariables = len(r_cross_var)
    if nvariables > 0:
        # check that the variables pass the same selection as the EDR3 stars
        # (this includes de colour cut)
        mask_var = []
        for item in r_cross_var['dr3_source_id']:
            if item in r_edr3['source_id']:
                mask_var.append(True)
            else:
                mask_var.append(False)
        r_cross_var = r_cross_var[mask_var]
        nvariables = len(r_cross_var)
        if verbose:
            r_cross_var.pprint(max_width=1000)
    else:
        r_cross_var = None

    return r_cross_var, nvariables
