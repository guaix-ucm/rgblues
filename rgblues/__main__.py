# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
RGB predictions of Gaia EDR3 stars

This code is hosted at https://github.com/guaix-ucm/rgblues
Authors: Nicolás Cardiel <cardiel@ucm.es>
         Sergio Pascual <sergiopr@fis.ucm.es>
         Rafael González <rafael08@ucm.es>

Usage example:
$ rgblues 56.66 24.10 1.0 12
"""

import argparse
import sys

from astropy.io import fits
import pooch

from .core.step1 import step1
from .core.step2 import step2
from .core.step3 import step3
from .core.step4 import step4
from .core.step5 import step5
from .core.step6 import step6
from .core.step7 import step7
from .gui.step8 import step8
from .version import version

MAX_SEARCH_RADIUS = 30  # degrees
EDR3_SOURCE_ID_15M_ALLSKY = 'edr3_source_id_15M_allsky.fits'
RGB_FROM_GAIA_ALLSKY = 'rgb_from_gaia_allsky.fits'


def right_ascension(ra_str):
    ra = float(ra_str)
    if not (0 <= ra <= 360):
        print('Right ascension must be 0 <= ra <= 360 degree')
        msg = f'Right ascension {ra} degree out of range'
        print(msg)
        raise ValueError(msg)
    return ra


def declination(dec_str):
    dec = float(dec_str)
    if not (-90 <= dec <= 90):
        print('Declination must be -90 <= dec <= 360 degree')
        msg = f'Declination {dec} degree out of range'
        print(msg)
        raise ValueError(msg)
    return dec


def search_radius(r_str):
    r = float(r_str)
    if not (0 < r < MAX_SEARCH_RADIUS):
        print(f'Search radius must be 0 < r <= {MAX_SEARCH_RADIUS} degree')
        msg = f'Search radius {r} degree out of range'
        print(msg)
        raise ValueError(msg)
    return r


def exec_rgblues(args):
    """Callable function that executes all the required steps

    Parameters
    ----------
    args: argparse instance
        Argparse instance containing the command-line parameters

    """
    print(f'\n        Welcome to rgblues version {version}')
    print(f'        ==============================\n')

    # check whether the auxiliary FITS binary table exists
    if args.debug:
        auxbintable = RGB_FROM_GAIA_ALLSKY
        auxhash = "md5:5e42a58471c780ef622d7bd620af3ea2"
    else:
        auxbintable = EDR3_SOURCE_ID_15M_ALLSKY
        auxhash = "md5:927f0dc8e74562268ee98bb13f6f88e3"

    fauxbin = pooch.retrieve(
        f"http://nartex.fis.ucm.es/~ncl/rgbphot/gaia/{auxbintable}",
        known_hash=auxhash
    )

    # read the previous file
    try:
        with fits.open(fauxbin) as hdul_table:
            edr3_source_id_15M_allsky = hdul_table[1].data.source_id
            if args.debug:
                table15M = hdul_table[1].data.copy()
            else:
                table15M = None
    except FileNotFoundError:
        raise SystemExit(f'ERROR: unexpected problem while reading {EDR3_SOURCE_ID_15M_ALLSKY}')

    # ---
    # step 1: EDR3 query
    r_edr3, nstars, nstars_colorcut = step1(
        args.ra_center,
        args.dec_center,
        args.search_radius,
        args.g_limit,
        args.verbose
    )

    # ---
    # step 2: retrieve additional parameters from the StarHorse sample
    if args.starhorse_block > 0:
        r_edr3 = step2(
            r_edr3,
            args.starhorse_block,
            args.verbose,
            args.debug
        )
    else:
        print('<STEP2> Retrieving StarHorse data from Gaia@AIP... (skipped!)')

    # ---
    # step 3: intersection with 15M star sample
    intersection = step3(
        r_edr3,
        edr3_source_id_15M_allsky,
        args.debug
    )

    # ---
    # step 4: DR2 query to identify variable stars
    r_dr2, nvariables, mask_var = step4(
        args.ra_center,
        args.dec_center,
        args.search_radius,
        args.g_limit,
        args.verbose
    )
    # ---

    # step 5: cross-match between DR2 and EDR3 to identify the variable stars
    if nvariables > 0:
        r_cross_var, nvariables = step5(
            r_dr2,
            mask_var,
            r_edr3,
            args.verbose
        )
    else:
        r_cross_var = None
        nvariables = 0
    print(f'        --> {nvariables} variable(s) in selected EDR3 star sample')

    # ---
    # step 6: compute RGB magnitudes
    r_edr3 = step6(
        r_edr3,
        args.verbose
    )

    # ---
    # step 7: generate output CSV files
    step7(
        r_edr3,
        args.basename,
        args.starhorse_block,
        nvariables,
        r_cross_var,
        intersection,
        table15M,
        args.verbose,
        args.debug
    )

    # ---
    # step 8: generate PDF chart
    if args.noplot:
        sys.stdout.write('<STEP8> No PDF plot generated (skipped!)\n')
        sys.stdout.flush()
    else:
        step8(
            r_edr3,
            args.ra_center,
            args.dec_center,
            args.search_radius,
            args.symbsize,
            args.brightlimit,
            args.nonumbers,
            args.nocolor,
            args.basename,
            nstars_colorcut,
            nvariables,
            r_cross_var,
            intersection,
            version,
            args.verbose
        )
    print('End of program')


def main():
    """Main function to parse input arguments"""
    parser = argparse.ArgumentParser(description=f"RGB predictions for Gaia EDR3 stars (version {version})")
    parser.add_argument("ra_center", help="right Ascension (decimal degrees)", type=right_ascension)
    parser.add_argument("dec_center", help="declination (decimal degrees)", type=declination)
    parser.add_argument("search_radius", help="search radius (decimal degrees)", type=search_radius)
    parser.add_argument("g_limit", help="limiting Gaia G magnitude", type=float)
    parser.add_argument("--basename", help="file basename for output files", type=str, default="rgblues")
    parser.add_argument("--brightlimit",
                        help="stars brighter than this Gaia G limit are displayed with star symbols (default=8.0)",
                        type=float, default=8.0)
    parser.add_argument("--symbsize", help="multiplying factor for symbol size (default=1.0)",
                        type=float, default=1.0)
    parser.add_argument("--nonumbers", help="do not display star numbers in PDF chart", action="store_true")
    parser.add_argument("--noplot", help="skip PDF chart generation", action="store_true")
    parser.add_argument("--nocolor", help="do not use colors in PDF chart", action="store_true")
    parser.add_argument("--starhorse_block", help="number of stars/query (default=0, no query)",
                        default=0, type=int)
    parser.add_argument("--verbose", help="increase program verbosity", action="store_true")
    parser.add_argument("--debug", help="debug flag", action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_usage()
        raise SystemExit()

    exec_rgblues(args)


if __name__ == "__main__":

    main()
