# RGBfromGaiaEDR3

This Python script predicts RGB magnitudes from *Gaia* EDR3 
photometric data. The code performs a cone search defined by coordinates 
right ascension and declination on the sky and a search radius. The 
predictions make use of the polynomial transformations given by Eqs. (2)-(5)
in [Cardiel et al. (2021; hereafter C21)](#1)

The cone search is performed making use of the Astroquery coordinated 
package of astropy. 

## Downloading and executing the code

In order to keep everything as simple as possible, this code does not need 
any particular installation. Simply download the script and execute it from 
the command line:

```buildoutcfg
$ git clone https://github.com/nicocardiel/RGBfromGaiaEDR3.git

$ cd RGBfromGaiaEDR3
$ python RGBfromGaiaEDR3.py 56.66 24.10 1.0 12
```
The last instruction executes the `RGBfromGaiaEDR3.py`script providing the 
four positional arguments: right ascension, declination, search radius and 
limiting *Gaia* G magnitude. *Note that the coordinates and search radius 
must be given in decimal degrees*.

The first time you execute the code, the auxiliary file 
`edr3_source_id_15M_allsky.fits`, containing the `source_id`of the *Gaia* 
EDR3 stars belonging to the ~15 million star sample of C21, is automatically 
downloaded. 

The script executes the following steps:

- Step 1: cone search in *Gaia* EDR3, gathering the following parameters: 
  `source_id`, `ra`, `dec`, `phot_g_mean_mag`, `phot_bp_mean_mag` and 
  `phot_rp_mean_mag`.
  
- Step 2: cross-matching of the previous EDR3 sample with the list of ~15 
  million stars from C21. This step determines the 
  subsample of EDR3 stars for which the RGB photometric calibration is 
  reliable.
  
- Step 3: cone search in *Gaia* DR2. This additional step is performed in 
  order to retrieve the `phot_variable_flag` parameter indicating whether 
  the star was flagged as variable in DR2. Note that this flag is not 
  available in EDR3.
  
- Step 4: cross-matching between DR2 and EDR3 to identify the variable 
  stars in EDR3. This step is required because it is not guaranteed that 
  the same astronomical source will always have the same source identifier 
  in the different Gaia Data Releases.
  
- Step 5: computation of the RGB magnitudes using the polynomial 
  transformations given in Eqs. (2)-(5) of C21.

- Step 6: generation of the output files. Three files (in CSV format) are 
  generated: 

    - `rgbsearch_15m.csv`: stars belonging to the ~15 million star sample 
      of C21 (with reliable RGB magnitude estimates)
      
    - `rgbsearch_var.csv`: stars flagged as variable in DR2.
    
    - `rgbsearch_edr3.csv`: remaining stars in EDR3. The RGB magnitude 
      estimates of these objects can be potentially biased due to 
      systematic effects introduced by interstellar extinction and 
      non-solar metallicity. This file will typically contain many more stars 
      than the `rgbsearch_15m.csv` selection.
      
  The three CSV files provide the same 10 parameters: 
  
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
    
- Step 7: generation of a finding chart plot (in PDF format): `rgbsearch.
  pdf`. The execution of the previous example generates a cone search around 
  the Pleiades star cluster:
  ![Pleiades chart](pleiades.png)
  The stars in this plot are color coded based on the *Gaia* G_BP - G_RP 
  colour. A red circle has been overplotted on the stars belonging to 
  the ~15 million star sample of C21, and a blue square on the variable 
  objects in DR2. Stars brighter than a pre-defined threshold are displayed 
  with big star symbols. 

Note that the four output files share the same root name `rgbsearch`. This 
can be easily modified using the optional argument `--basename 
<newbasename>` in the command line.

### Additional help

Some auxiliary optional arguments are also available. See description 
invoking the script help:

```buildoutcfg
$ python RGBfromGaiaEDR3.py --help

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
  --noplot              skip PDF chart generation
  --nocolor             do not use colors in PDF chart
  --verbose             increase program verbosity
```

## Citation
If you find this Python script useful, please cite:

<a id="1">Cardiel et al. (2021)</a>, MNRAS, in preparation