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
#### Input files

This is a .dat file where the header is either "F444W F200W-F444W Separation" (for NIRCam) or "F1500W F1500W-F2100W Separation" (for MIRI) 
* see input_miri_example.dat/input_nircam_example.dat for the format
* separation is in arcseconds

### Running the dropout test
* This is optional as I have already done this (see the "keep" photometry file in JADES and SMILES folders)
* If you want to run this for yourself:

```
python dropout_test.py 'fits file directory' Instrument output_file_name.dat
```

You can download the relavant fits files at https://archive.stsci.edu/hlsp/jades (JADES goods-s field: F200W, F356W, F444W) and https://archive.stsci.edu/hlsp/smiles (SMILES goods-s field: F560W, F1500W, F2100W)

## Author:

Aiza Kenzhebekova

a.kenzhebekova@sms.ed.ac.uk
