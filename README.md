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

## Running the dropout test
* This is optional as I have already done this (see the "keep" photometry file in JADES and SMILES folders)
* If you want to run this for yourself:

```
python dropout_test.py 'fits file directory' Instrument output_file_name.dat
```
You can download the relavant fits files at https://archive.stsci.edu/hlsp/jades (JADES goods-s field: F356W, F200W, F444W) and 

## Author:

Aiza Kenzhebekova
a.kenzhebekova@sms.ed.ac.uk
