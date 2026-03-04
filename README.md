# SIDERIS
Statistical Identification of Directly-imaged Exoplanet candidates, Reviewed to Inform Selection

## Description
SIDERIS calculates the probability of finding a background object with specified colour and magnitude within a separation from a target star.

## Getting Started

### Running the calculator
To run the background probabiliy calculator:
```
python background_object_prob_calc.py input_miri_example.dat
```

or if you would like to run TRILEGAL (only recomended if you have blue colour candidates):
```
python background_object_prob_calc.py input_miri_example.dat --trilegal RA Dec FOV
```
* Sometimes, TRILEGAL fails to run. If this happens, please try again at a later time, or run without TRILEGAL.
* RA and Dec in decimal degrees
* FOV (field of view) in sq. arcsec
  
#### Input files

This is a .dat file where the header is either "F444W F200W-F444W Separation" (for NIRCam) or "F1500W F1500W-F2100W Separation" (for MIRI) 
* see input_miri_example.dat/input_nircam_example.dat for the format
* separation is in arcseconds

### Running the dropout test
* This is optional as I have already done this (see the "keep" photometry file in JADES and SMILES folders)
* If you want to run this:
```
python dropout_test.py 'fits file directory' Instrument output_file_name.dat
```

You can download the relavant fits files at:
* https://archive.stsci.edu/hlsp/jades (JADES goods-s field: F200W, F356W, F444W)
* https://archive.stsci.edu/hlsp/smiles (SMILES goods-s field: F560W, F1500W, F2100W)

## References and Acknowledgements
Alberts S., et al., 2024, ApJ, 976, 224 doi:[10.3847/1538-4357/ad7396](https://iopscience.iop.org/article/10.3847/1538-4357/ad7396)

Astropy Collaboration, Price-Whelan, A. M., Lim, P. L., et al. 2022, ApJ, 935, 167 doi:[10.3847/1538-4357/ac7c74](https://iopscience.iop.org/article/10.3847/1538-4357/ac7c74)

Eisenstein D. J., et al., 2023, arXiv e-prints, p. arXiv:2306.02465 doi:[10.3847/1538-4365/ae3163](https://iopscience.iop.org/article/10.3847/1538-4365/ae3163)

Linder E. F., Mordasini C., Mollière P., Marleau G.-D., Malik M., Quanz S. P., Meyer M. R., 2019, A&A, 623, A85 doi:[10.1051/0004-6361/201833873](https://www.aanda.org/articles/aa/full_html/2019/03/aa33873-18/aa33873-18.html)

Rieke G. H., Alberts S., Shivaei I., Lyu J., Willmer C. N. A., Pérez-González P., Williams C. C., 2024, ApJ, 975, 83 doi:[10.3847/1538-4357/ad6cd2](https://iopscience.iop.org/article/10.3847/1538-4357/ad6cd2)

This work is based [in part] on observations made with the NASA/ESA/CSA James Webb Space Telescope. The data were obtained from the Mikulski Archive for Space Telescopes at the Space Telescope Science Institute, which is operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract NAS 5-03127 for JWST.

## Author:

Aiza Kenzhebekova

Email: a.kenzhebekova@sms.ed.ac.uk
