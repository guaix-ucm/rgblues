# -*- coding: utf-8 -*-
#
# Copyright 2021 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Test generation of pleiades cone search used in the paper
"""

import pooch

from rgblues.__main__ import exec_rgblues
from rgblues.core.step7 import OUTTYPES


class Args:
    ra_center = 56.66
    dec_center = 24.10
    search_radius = 1.0
    g_limit = 12.0
    basename = 'rgblues'
    brightlimit = 8.0
    symbsize = 1.0
    nonumbers = False
    noplot = True            # do not generate PDF file
    nocolor = False
    starhorse_block = 0
    verbose = False
    debug = False


auxhash = {
    'edr3': "md5:5b443fbc2863ff8c9c10bfec791acfb2",
    '15m': "md5:34a416e7e25235d50ba609d2e94e8a49",
    'var': "md5:81fde79fabb87cf0941c875b6e23c03d"
}

fref = dict()
print('Using reference files:')
for ftype in OUTTYPES:
    fname = f'{Args.basename}_{ftype}.csv'
    ftmp = pooch.retrieve(
        f"http://nartex.fis.ucm.es/~ncl/rgbphot/gaia/{fname}",
        known_hash=auxhash[ftype]
    )
    fref[ftype] = ftmp
    print(ftmp)


def test_pleiades_simple():
    exec_rgblues(Args)
    for ftype in OUTTYPES:
        fname = f'{Args.basename}_{ftype}.csv'
        print(f'Checking {fname}')
        with open(fname, 'r') as f:
            file1 = f.read()
        with open(fref[ftype], 'r') as f:
            file2 = f.read()
        assert file1 == file2
