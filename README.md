# rgblues

## Introduction

This Python package `rgblues` predicts RGB magnitudes from *Gaia* EDR3 
photometric data. These magnitudes are given in the standard system defined by
[Cardiel et al. (2021a)](#1).

The code performs a cone search defined by coordinates 
right ascension and declination on the sky and a search radius. The 
predictions make use of the polynomial transformations given by Eqs. (2)-(5)
in [Cardiel et al. (2021b; hereafter C21)](#2)

The cone search is performed making use of the Astroquery coordinated 
package of astropy. 

You need to have a live connection to the Internet for 
the script to work!

**Note**: You might also be interested in the Python package 
[rgbloom](https://pypi.org/project/rgbloom/), which is based on the work
published by [Carrasco et al. (2023)](#3).

## Installing the code

In order to keep your current Python installation clean, it is highly 
recommended to install a python 3 *virtual environment* first.

### Creating and activating the python virtual environment

```bash
$ python3 -m venv venv_rgb
$ . venv_rgb/bin/activate
(venv_rgb) $
```

### Installing the package

We recommend installing the latest stable version, which is available via the
[PyPI repository](https://pypi.org/project/rgblues/):

```bash
(venv_rgb) $ pip install rgblues
```

The latest development version is available through [GitHub](https://github.com/nicocardiel/rgblues):

```bash
(venv_rgb) $ pip install git+https://github.com/guaix-ucm/rgblues.git@main#egg=rgblues
```

## Executing the program
Just execute it from the command line:

```buildoutcfg
(venv_rgb) $ rgblues 56.66 24.10 1.0 12
```

The last instruction executes the program providing the 
four positional arguments: right ascension, declination, search radius and 
limiting *Gaia* G magnitude. *Note that the coordinates and search radius 
must be given in decimal degrees*.

The first time you execute the code, the auxiliary file
`edr3_source_id_15M_allsky.fits` (size 129 Mb), containing the `source_id`of
the *Gaia* EDR3 stars belonging to the ~15 million star sample of C21, is
automatically downloaded to a cache directory (you do not have to worry
about its location). 

The execution of this example should led to the following output in the terminal 
(except for the absolute path where the auxiliary downloaded file is stored):

```
        Welcome to rgblues version 1.2
        ==============================

Downloading data from 'http://nartex.fis.ucm.es/~ncl/rgbphot/gaia/edr3_source_id_15M_allsky.fits' to file '/Users/cardiel/Library/Caches/pooch/bf659d7a02408a54b3fa6b62fb3b051e-edr3_source_id_15M_allsky.fits'.
<STEP1> Starting cone search in Gaia EDR3... (please wait)
  INFO: Query finished. [astroquery.utils.tap.core]
        --> 310 stars found
        --> 19 stars outside -0.5 < G_BP-G_RP < 2.0
<STEP2> Retrieving StarHorse data from Gaia@AIP... (skipped!)
<STEP3> Cross-matching EDR3 with 15M subsample... (please wait)
        --> 55 stars in common with 15M sample
<STEP4> Looking for variable stars in Gaia DR2... (please wait)
  INFO: Query finished. [astroquery.utils.tap.core]
        --> 310 stars in DR2, (2 initial variables)
<STEP5> Cross-matching variables in DR2 with stars in EDR3... (please wait)
  INFO: Query finished. [astroquery.utils.tap.core]
        --> 2 variable(s) in selected EDR3 star sample
<STEP6> Computing RGB magnitudes...OK
<STEP7> Saving output CSV files...OK
<STEP8> Generating PDF plot...OK
End of program
```

The script executes the following steps:

- Step 1: cone search in *Gaia* EDR3, gathering the following parameters: 
  `source_id`, `ra`, `dec`, `phot_g_mean_mag`, `phot_bp_mean_mag` and 
  `phot_rp_mean_mag`.

- Step 2: cone search in StarHorse to retrieve interstellar extinction,
  metallicity and distance, among other parameters. This step is optional and
  only executed when `--starhorse_block <number>` is employed (in this case
  `<number>` is an integer number indicating the number of stars whose
  parameters are retrieved in each single query to Gaia@AIP; a typical useful 
  value is 100).
  
- Step 3: cross-matching of the previous EDR3 sample with the list of ~15 
  million stars from C21. This step determines the 
  subsample of EDR3 stars for which the RGB photometric calibration is 
  reliable.
  
- Step 4: cone search in *Gaia* DR2. This additional step is performed in 
  order to retrieve the `phot_variable_flag` parameter indicating whether 
  the star was flagged as variable in DR2. Note that this flag is not 
  available in EDR3.
  
- Step 5: cross-matching between DR2 and EDR3 to identify the variable 
  stars in EDR3. This step is required because it is not guaranteed that 
  the same astronomical source will always have the same source identifier 
  in the different Gaia Data Releases.
  
- Step 6: computation of the RGB magnitudes using the polynomial 
  transformations given in Eqs. (2)-(5) of C21.

- Step 7: generation of the output files. Three files (in CSV format) are 
  generated: 

    - `rgblues_15m.csv`: stars belonging to the ~15 million star sample 
      of C21 (with reliable RGB magnitude estimates).

    - `rgblues_var.csv`: objects flagged as variable in DR2.
    
    - `rgblues_edr3.csv`: remaining objects in EDR3. The RGB magnitude 
      estimates of these objects can be potentially biased due to 
      systematic effects introduced by interstellar extinction, or by 
      exhibiting non-solar metallicity, or a colour outside the *Gaia* -0.5 < 
      G_BP-G_RP < 2.0 interval. This file will typically contain more stars 
      than the `rgblues_15m.csv` selection.
      
  The three CSV files provide the same 11 columns:
  
    - `number`: consecutive number of the star in each CSV file
    - `source_id`: identification in EDR3
    - `ra`: right ascension (from EDR3)
    - `dec`: declination (from EDR3)
    - `b_rgb`: blue RGB magnitude estimate
    - `g_rgb`: green RGB magnitude estimate
    - `r_rgb`: red RGB magnitude estimate
    - `g_br_rgb`: pseudo-green RGB magnitude estimate, defined in C21 as 
      the arithmetic mean of the blue and red RGB magnitudes
    - `phot_g_mean_mag`: *Gaia* G magnitude (EDR3)
    - `phot_bp_mean_mag`: *Gaia* G_BP magnitude (EDR3)
    - `phot_rp_mean_mag`: *Gaia* G_RP magnitude (EDR3)

  The list of objects in those files is sorted by right ascension.

  When using `--starhorse_block <number>`, the files `rgblues_15m.csv` and
  `rgblues_edr3.csv` contain 3 additional
  columns providing parameters derived by [Anders et al. (2019)](#3):

    - `av50`: 50th percentile of the interstellar extinction 
    - `met50`: 50th percentile of the metallicity [M/H]
    - `dist50`: 50th percentile of the distance (kpc)

  These three values are set to 99.999 for those stars that do not belong to
  the StarHorse sample.

- Step 8: generation of a finding chart plot (in PDF format): `rgblues.pdf`. 
  The execution of the previous example generates a cone search around 
  the [Pleiades](https://en.wikipedia.org/wiki/Pleiades) star cluster:
  ![Pleiades plot](http://nartex.hst.ucm.es/~ncl/rgbphot/gaia/pleiades_v5.png)
  The stars in this plot 
  (see [PDF file](http://nartex.hst.ucm.es/~ncl/rgbphot/gaia/pleiades_v5.pdf))
  are color coded based on the *Gaia* G_BP - G_RP 
  colour. A red circle has been overplotted on the stars belonging to 
  the ~15 million star sample of C21, a blue square on the variable 
  objects in DR2, and a grey diamond on EDR3 stars outside the *Gaia* 
  -0.5 < G_BP - G_RP < 2.0 colour interval. 
  Stars brighter than a pre-defined threshold are displayed 
  with big star symbols. To facilitate the identification of each star, the
  consecutive star number in the three files (`rgblues_15m.csv`,
  `rgblues_edr3.csv` and `rgblues_var.csv`) is also displayed (in red,
  black and blue, respectively). These numbers are not displayed when using the
  parameter `--nonumbers` in the command line.

Note that the four output archives (1 PDF and 3 CSV files) share the same root
name `rgblues`. This can be easily modified using the optional argument
`--basename <newbasename>` in the command line.

### Additional help

Some auxiliary optional arguments are also available. See description 
invoking the script help:

```buildoutcfg
$ rgblues --help

...
...

positional arguments:
  ra_center             right Ascension (decimal degrees)
  dec_center            declination (decimal degrees)
  search_radius         search radius (decimal degrees)
  g_limit               limiting Gaia G magnitude

optional arguments:
  -h, --help            show this help message and exit
  --basename BASENAME   file basename for output files
  --brightlimit BRIGHTLIMIT
                        stars brighter than this Gaia G limit are displayed 
                        with star symbols (default=8.0)
  --symbsize SYMBSIZE   multiplying factor for symbol size (default=1.0)
  --nonumbers           do not display star numbers in PDF chart
  --noplot              skip PDF chart generation
  --nocolor             do not use colors in PDF chart
  --starhorse_block STARHORSE_BLOCK
                        number of stars/query (default=0, no query)
  --verbose             increase program verbosity
  --debug               debug flag
```

## Citation
If you find this Python package useful, 
please cite [Cardiel et al. (2021a)](#1)
(to quote the use of the standard RGB system)
and [Cardiel et al. (2021b)](#2) (where the transformation between the *Gaia*
photometry and the RGB magnitudes is derived).

## Related information

You can visit the [RGB Photometry](https://guaix.ucm.es/rgbphot) web page at
the Universidad Complutense de Madrid.

## Bibliography

<a id="3">Anders et al. (2019)</a>, 
https://ui.adsabs.harvard.edu/abs/2019A%26A...628A..94A/abstract

<a id="1">Cardiel et al. (2021a)</a>, 
MNRAS, https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.3730C/abstract

<a id="2">Cardiel et al. (2021b)</a>, 
MNRAS, https://ui.adsabs.harvard.edu/abs/2021MNRAS.507..318C/abstract

<a id="3">Carrasco et al. (2023)</a>, Remote Sensing, https://www.mdpi.com/2072-4292/15/7/1767
