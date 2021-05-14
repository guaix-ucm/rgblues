# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
RGB predictions of Gaia EDR3 stars

This code is hosted at https://github.com/nicocardiel/RGBfromGaiaEDR3
Maintainer: Nicolás Cardiel <cardiel@ucm.es>

Example:
    python RGBfromGaiaEDR3.py 56.66398954459 24.10299456629 1 12
"""

import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Column
from astropy.wcs import WCS
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 17.0
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 2
import numpy as np
from numpy.polynomial import Polynomial
import os
import sys
import urllib

MAX_SEARCH_RADIUS = 30  # degrees
EDR3_SOURCE_ID_15M_ALLSKY = 'edr3_source_id_15M_allsky.fits'
VERSION = 1.0

def main():

    parser = argparse.ArgumentParser(description="RGB predictions for Gaia EDR3 stars")
    parser.add_argument("ra_center", help="Right Ascension (decimal degrees)", type=float)
    parser.add_argument("dec_center", help="Declination (decimal degrees)", type=float)
    parser.add_argument("search_radius", help="Search radius (decimal degrees)", type=float)
    parser.add_argument("g_limit", help="Limiting Gaia G magnitude", type=float)
    parser.add_argument("--basename", help="File basename for output files", type=str, default="rgbchart")
    parser.add_argument("--brightlimit",
                        help="Stars brighter that this Gaia G limit are displayed with star symbols",
                        type=float, default=8.0)
    parser.add_argument("--nocolor", help="Do not use colors in PDF chart", action="store_true")
    parser.add_argument("--verbose", help="Increase program verbosity", action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_usage()
        raise SystemExit()

    if args.ra_center < 0 or args.ra_center > 360:
        raise SystemExit('ERROR: right ascension out of valid range')
    if args.dec_center < -90 or args.dec_center > 90:
        raise SystemExit('ERROR: declination out of valid range')
    if args.search_radius < 0:
        raise SystemExit('ERROR: search radius must be > 0 degrees')
    if args.search_radius > MAX_SEARCH_RADIUS:
        raise SystemExit(f'ERROR: search radius must be <= {MAX_SEARCH_RADIUS} degrees')

    # check whether the FITS file containing the source_id of the 15M star sample exists
    if os.path.isfile(EDR3_SOURCE_ID_15M_ALLSKY):
        pass
    else:
        urldir = f'http://nartex.fis.ucm.es/~ncl/rgbphot/gaia/{EDR3_SOURCE_ID_15M_ALLSKY}'
        sys.stdout.write(f'Downloading {urldir}... (please wait)')
        sys.stdout.flush()
        urllib.request.urlretrieve(urldir, EDR3_SOURCE_ID_15M_ALLSKY)
        print(' ...OK!')

    # read the previous file
    try:
        with fits.open(EDR3_SOURCE_ID_15M_ALLSKY) as hdul_table:
            edr3_source_id_15M_allsky = hdul_table[1].data.source_id
    except FileNotFoundError:
        raise SystemExit(f'ERROR: unexpected problem while reading {EDR3_SOURCE_ID_15M_ALLSKY}')

    # define WCS
    naxis1 = 1024
    naxis2 = naxis1
    pixscale = 2 * args.search_radius / naxis1

    wcs_image = WCS(naxis=2)
    wcs_image.wcs.crpix = [naxis1 / 2, naxis2 / 2]
    wcs_image.wcs.crval = [args.ra_center, args.dec_center]
    wcs_image.wcs.cunit = ["deg", "deg"]
    wcs_image.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_image.wcs.cdelt = [-pixscale, pixscale]
    wcs_image.array_shape = [naxis1, naxis2]
    if args.verbose:
        print(wcs_image)

    # EDR3 query
    query = f"""
    SELECT source_id, ra, dec,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag

    FROM gaiaedr3.gaia_source
    WHERE 1=CONTAINS(
      POINT('ICRS', {args.ra_center}, {args.dec_center}), 
      CIRCLE('ICRS',ra, dec, {args.search_radius}))
    AND phot_g_mean_mag IS NOT NULL 
    AND phot_bp_mean_mag IS NOT NULL 
    AND phot_rp_mean_mag IS NOT NULL
    AND phot_bp_mean_mag - phot_rp_mean_mag > -0.5
    AND phot_bp_mean_mag - phot_rp_mean_mag < 2.0
    AND phot_g_mean_mag < {args.g_limit}
    
    ORDER BY ra
    """
    sys.stdout.write('Starting Gaia EDR3 query... (please wait) ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_edr3 = job.get_results()
    nstars = len(r_edr3)
    print(f'--> {nstars} stars found\n')
    if nstars == 0:
        raise SystemExit('ERROR: no stars found. Change search parameters!')
    if args.verbose:
        r_edr3.pprint(max_width=1000)

    # intersection with 15M star sample
    sys.stdout.write('Cross-matching EDR3 with 15M subsample... (please wait) ')
    sys.stdout.flush()
    set1 = set(np.array(r_edr3['source_id']))
    set2 = set(edr3_source_id_15M_allsky)
    intersection = set2.intersection(set1)
    print('')
    if args.verbose:
        print(len(set1), len(set2), len(intersection))

    # DR2 query to identify variable stars
    query = f"""
    SELECT source_id, ra, dec, phot_g_mean_mag, phot_variable_flag

    FROM gaiadr2.gaia_source
    WHERE  1=CONTAINS(
      POINT('ICRS', {args.ra_center}, {args.dec_center}), 
      CIRCLE('ICRS',ra, dec, {args.search_radius}))
    AND phot_g_mean_mag < {args.g_limit}
    """
    sys.stdout.write('Starting Gaia DR2 query... (please wait) ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_dr2 = job.get_results()
    print('')
    nstars_dr2 = len(r_dr2)
    if nstars_dr2 == 0:
        nvariables = 0
        mask = None
    else:
        print(f'Number of DR2 objects: {nstars_dr2}')
        if isinstance(r_dr2['phot_variable_flag'][0], bytes):
            mask = r_dr2['phot_variable_flag'] == b'VARIABLE'
        elif isinstance(r_dr2['phot_variable_flag'][0], str):
            mask = r_dr2['phot_variable_flag'] == 'VARIABLE'
        else:
            raise SystemExit('Unexpected type of data in column phot_variable_flag')
        nvariables = sum(mask)
        print(f'Number of variable stars: {nvariables}')
    if nvariables > 0:
        if args.verbose:
            r_dr2[mask].pprint(max_width=1000)

    # cross-match between DR2 and EDR3 to identify the variable stars
    dumstr = '('
    if nvariables > 0:
        # generate sequence of source_id of variable stars
        for i, item in enumerate(r_dr2[mask]['source_id']):
            if i > 0:
                dumstr += ', '
            dumstr += f'{item}'
        dumstr += ')'
        # cross-match
        query = f"""
        SELECT *
        FROM gaiaedr3.dr2_neighbourhood
        WHERE dr2_source_id IN {dumstr}
        ORDER BY angular_distance
        """
        sys.stdout.write('Starting query in gaiaedr3.dr2_neighbourhood... (please wait) ')
        sys.stdout.flush()
        job = Gaia.launch_job_async(query)
        r_cross = job.get_results()
        print('')
        if args.verbose:
            r_cross.pprint(max_width=1000)

        # compute G_BP - G_RP colour and predict RGB magnitudes
        r_edr3.add_column(
            Column(r_edr3['phot_bp_mean_mag'] - r_edr3['phot_rp_mean_mag'],
                   name='bp_rp', unit=u.mag)
        )

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
        if args.verbose:
            r_edr3.pprint(max_width=1000)

        # ToDo: segregate output (stars in 15M sample, variables, remaining in EDR3
        # save result in CSV file
        ascii.write(r_edr3, f'{args.basename}.csv', overwrite=True,
                    include_names=['source_id', 'ra', 'dec', 'b_rgb', 'g_rgb', 'r_rgb', 'g_br_rgb',
                                   'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag'],
                    format='csv', delimiter=',',
                    formats={'source_id': '19d',
                             'ra': '14.9f',
                             'dec': '14.9f',
                             'b_rgb': '6.2f',
                             'g_rgb': '6.2f',
                             'r_rgb': '6.2f',
                             'g_br_rgb': '6.2f',
                             'phot_g_mean_mag': '8.4f',
                             'phot_bp_mean_mag': '8.4f',
                             'phot_rp_mean_mag': '8.4f'
                             })

        # ---

        # generate plot
        r_edr3.sort('phot_g_mean_mag')
        if args.verbose:
            r_edr3.pprint(max_width=1000)

        symbol_size = (50 / np.array(r_edr3['phot_g_mean_mag'])) ** 2.5
        ra_array = np.array(r_edr3['ra'])
        dec_array = np.array(r_edr3['dec'])

        c = SkyCoord(ra=ra_array * u.degree, dec=dec_array * u.degree, frame='icrs')
        x_pix, y_pix = wcs_image.world_to_pixel(c)

        fig = plt.figure(figsize=(11, 10))
        ax = plt.subplot(projection=wcs_image)
        iok = r_edr3['phot_g_mean_mag'] < args.brightlimit
        if args.nocolor:
            sc = ax.scatter(x_pix[iok], y_pix[iok], marker='*', color='grey',
                       edgecolors='black', linewidth=0.2, s=symbol_size[iok])
            ax.scatter(x_pix[~iok], y_pix[~iok], marker='.', color='grey',
                       edgecolors='black', linewidth=0.2, s=symbol_size[~iok])
        else:
            cmap = plt.cm.get_cmap('jet')
            sc= ax.scatter(x_pix[iok], y_pix[iok], marker='*',
                       edgecolors='black', linewidth=0.2, s=symbol_size[iok],
                       cmap=cmap, c=r_edr3[iok]['bp_rp'])
            ax.scatter(x_pix[~iok], y_pix[~iok], marker='.',
                       edgecolors='black', linewidth=0.2, s=symbol_size[~iok],
                       cmap=cmap, c=r_edr3[~iok]['bp_rp'])

        sorter = np.argsort(r_edr3['source_id'])
        iok = np.array(sorter[np.searchsorted(r_edr3['source_id'], r_dr2[mask]['source_id'], sorter=sorter)])
        ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='s', facecolors='none', edgecolors='blue', linewidth=0.5)

        sorter = np.argsort(r_edr3['source_id'])
        iok = np.array(sorter[np.searchsorted(r_edr3['source_id'], np.array(list(intersection)), sorter=sorter)])
        ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='o', facecolors='none', edgecolors='red', linewidth=0.5)

        ax.set_xlabel('ra')
        ax.set_ylabel('dec')

        ax.set_aspect('equal')

        if not args.nocolor:
            cbaxes = fig.add_axes([0.20, 0.84, 0.15, 0.02])
            cbar = plt.colorbar(sc, cax=cbaxes, orientation='horizontal', format='%3.1f')
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label(label=r'$G_{\rm BP}-G_{\rm RP}$', size=12)

        overlay = ax.get_coords_overlay('icrs')
        overlay.grid(color='black', ls='dotted')

        ax.text(0.99, 0.01, f'RGBfromGaiaEDR3, version {VERSION}', fontsize=12,
                horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

        plt.savefig(f'{args.basename}.pdf')
        plt.close(fig)


if __name__ == "__main__":

    main()