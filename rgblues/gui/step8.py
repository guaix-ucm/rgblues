# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Generate output PDF plot
"""

import sys

import matplotlib.pyplot as plt
import matplotlib.style
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np

from .style import mpl_style
from rgblues.core.step7 import OUTTYPES
OUTTYPES_COLOR = {'edr3': 'black', '15m': 'red', 'var': 'blue'}


def step8(r_edr3, ra_center, dec_center, search_radius, symbsize, brightlimit,
          nonumbers, nocolor, basename, nstars_colorcut, nvariables,
          r_cross_var, intersection, version, verbose):
    """Perform EDR3 query

    Parameters
    ----------
    r_edr3 : astropy Table
        Table containing the EDR3 query.
    ra_center : float
        Right ascension (decimal degree) corresponding to the center
        of the field of view.
    dec_center : float
        Declination (decimal degree) corresponding to the center
        of the field of view.
    search_radius : float
        Radius (decimal degrees) of the field of view.
    symbsize : float
        Multiplying factor for symbol size.
    brightlimit : float
        Stars brighter than this Gaia G limit are displayed with star
        symbols.
    nonumbers : bool
        If True, do not display star numbers in PDF chart.
    nocolor : bool
        If True, do not use colors in PDF chart.
    basename : str
        Base name for output files.
    nstars_colorcut : int
        Number of stars outside the colour cut
        -0.5 < G_BP - G_RP < 2.0 mag.
    nvariables : int
        Final number of variable stars in EDR3 query.
    r_cross_var : astropy Table
        Table with the cross-match of the DR2 and EDR3 queries.
    intersection : numpy array
        Array with the source_id values of the intersection between
        the EDR3 query and the 15M star sample.
    version : str
        Version number.
    verbose : bool
        If True, display additional information.

    """
    sys.stdout.write('<STEP8> Generating PDF plot...')
    sys.stdout.flush()

    # define WCS
    naxis1 = 1024
    naxis2 = naxis1
    pixscale = 2 * search_radius / naxis1

    wcs_image = WCS(naxis=2)
    wcs_image.wcs.crpix = [naxis1 / 2, naxis2 / 2]
    wcs_image.wcs.crval = [ra_center, dec_center]
    wcs_image.wcs.cunit = ["deg", "deg"]
    wcs_image.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_image.wcs.cdelt = [-pixscale, pixscale]
    wcs_image.array_shape = [naxis1, naxis2]
    if verbose:
        print(wcs_image)

    # generate plot
    r_edr3.sort('phot_g_mean_mag')
    if verbose:
        print('')
        r_edr3.pprint(max_width=1000)

    symbol_size = symbsize * (50 / np.array(r_edr3['phot_g_mean_mag'])) ** 2.5
    ra_array = np.array(r_edr3['ra'])
    dec_array = np.array(r_edr3['dec'])

    c = SkyCoord(ra=ra_array * u.degree, dec=dec_array * u.degree, frame='icrs')
    x_pix, y_pix = wcs_image.world_to_pixel(c)

    matplotlib.style.use(mpl_style)
    fig = plt.figure(figsize=(13, 10))
    ax = plt.subplot(projection=wcs_image)
    iok = r_edr3['phot_g_mean_mag'] < brightlimit
    if nocolor:
        sc = ax.scatter(x_pix[iok], y_pix[iok], marker='*', color='grey',
                        edgecolors='black', linewidth=0.2, s=symbol_size[iok])
        ax.scatter(x_pix[~iok], y_pix[~iok], marker='.', color='grey',
                   edgecolors='black', linewidth=0.2, s=symbol_size[~iok])
    else:
        cmap = plt.cm.get_cmap('jet')
        sc = ax.scatter(x_pix[iok], y_pix[iok], marker='*',
                        edgecolors='black', linewidth=0.2, s=symbol_size[iok],
                        cmap=cmap, c=r_edr3[iok]['bp_rp'], vmin=-0.5, vmax=2.0)
        ax.scatter(x_pix[~iok], y_pix[~iok], marker='.',
                   edgecolors='black', linewidth=0.2, s=symbol_size[~iok],
                   cmap=cmap, c=r_edr3[~iok]['bp_rp'], vmin=-0.5, vmax=2.0)

    # display numbers if requested
    if not nonumbers:
        for irow in range(len(r_edr3)):
            number_csv = r_edr3[irow]['number_csv']
            text = r_edr3[irow][f'number_{OUTTYPES[number_csv]}']
            ax.text(x_pix[irow], y_pix[irow], text,
                    color=OUTTYPES_COLOR[OUTTYPES[number_csv]], fontsize='5',
                    horizontalalignment='left', verticalalignment='bottom')

    # stars outside the -0.5 < G_BP - G_RP < 2.0 colour cut
    if nstars_colorcut > 0:
        mask_colour = np.logical_or((r_edr3['bp_rp'] <= -0.5), (r_edr3['bp_rp'] >= 2.0))
        iok = np.argwhere(mask_colour)
        ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='D', facecolors='none', edgecolors='grey', linewidth=0.5)

    # variable stars
    if nvariables > 0:
        sorter = np.argsort(r_edr3['source_id'])
        iok = np.array(sorter[np.searchsorted(r_edr3['source_id'], r_cross_var['dr3_source_id'], sorter=sorter)])
        ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='s', facecolors='none', edgecolors='blue', linewidth=0.5)

    # stars in 15M sample
    if len(intersection) > 0:
        sorter = np.argsort(r_edr3['source_id'])
        iok = np.array(sorter[np.searchsorted(r_edr3['source_id'], np.array(list(intersection)), sorter=sorter)])
        ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='o', facecolors='none',
                   edgecolors=OUTTYPES_COLOR['15m'], linewidth=0.5)

    ax.scatter(0.03, 0.96, s=240, marker='o', facecolors='white',
               edgecolors=OUTTYPES_COLOR['15m'], linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.96, 'star in 15M sample', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.scatter(0.03, 0.92, s=240, marker='s', facecolors='white',
               edgecolors=OUTTYPES_COLOR['var'], linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.92, 'variable in Gaia DR2', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.scatter(0.03, 0.88, s=240, marker='D', facecolors='white', edgecolors='grey', linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.88, 'outside colour range', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.set_xlabel('ra')
    ax.set_ylabel('dec')

    ax.set_aspect('equal')

    if not nocolor:
        cbaxes = fig.add_axes([0.683, 0.81, 0.15, 0.02])
        cbar = plt.colorbar(sc, cax=cbaxes, orientation='horizontal', format='%1.0f')
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(label=r'$G_{\rm BP}-G_{\rm RP}$', size=12, backgroundcolor='white')

    ax.text(0.98, 0.96, f'Field radius: {search_radius:.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.02, 0.06, r'$\alpha_{\rm center}$:', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.25, 0.06, f'{ra_center:.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.02, 0.02, r'$\delta_{\rm center}$:', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.25, 0.02, f'{dec_center:+.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.98, 0.02, f'rgblues, version {version}', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

    f = np.pi / 180
    xp = naxis1 / 2 + search_radius/pixscale * np.cos(np.arange(361)*f)
    yp = naxis2 / 2 + search_radius/pixscale * np.sin(np.arange(361)*f)
    ax.plot(xp, yp, '-', color='orange', linewidth=0.5, alpha=0.5)

    ax.set_xlim([-naxis1*0.12, naxis1*1.12])
    ax.set_ylim([-naxis2*0.05, naxis2*1.05])

    ax.set_axisbelow(True)
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls='dotted')

    plt.savefig(f'{basename}.pdf')
    plt.close(fig)
    if verbose:
        pass
    else:
        print('OK')
