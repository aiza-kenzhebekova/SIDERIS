""" 
Drop Out Test
MPhys Project 
Aiza Kenzhebekova 
02/03/2026
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from functions import dropout_test_functions as drptst 
from astropy.table import Table
import os
from pathlib import Path

def main():
    """ 
    This programme evaluates which sources drop out in the SMILES (MIRI) and JADES (NIRCam) galaxy surveys.
    --------------------------------------------------------------------------------------------------------
    Format: python drop_out_test.py "instrument" "outputfilename"

    Output: .dat file of the photometry and FWHM (in arcsec) of sources which do not drop out

    Directory: File path of the .fits files  
    Instrument: "NIRCam" or "MIRI"
    Suggested outfile name: "SMILES/JADES_photometry_retained.dat 
    """

    # check for correct inputs:
    if len(sys.argv) != 3:
        print("Usage: python drop_out_test.py <directory> <instrument> <outputfilename>")
        sys.exit(1)

    data_dir = sys.argv[1] # mine is '/Users/aizakenzhebekova/Downloads/MPhys_Project/
    instrument = sys.argv[2]
    outputfilename = sys.argv[3]

    BASE_DIR = Path(__file__).resolve().parent

    linder_file = os.path.join(BASE_DIR/'data_files'/ 'BEX_evol_mags_-2_MH_0.00.dat')

    # read in Linder+2019 model
    linder = Table.read(linder_file, format='ascii', names=('1:log(Age/yr)', '2:Mass/Mearth', '3:Radius/Rjupiter', 
                '4:Luminosity/Ljupiter', '5:Teff/K', '6:logg/cgs', '7:NACOJ', '8:NACOH', '9:NACOKs', '10:NACOLp', '11:NACOMp', '12:CousinsR', '13:CousinsI', '14:WISE1', '15:WISE2', '16:WISE3', 
                '17:WISE4', '18:F115W', '19:F150W', '20:F200W', '21:F277W', '22:F356W', '23:F444W', '24:F560W', '25:F770W', '26:F1000W', '27:F1280W', '28:F1500W', '29:F1800W', '30:F2100W', '31:F2550W',
                '32:VISIRB87', '33:VISIRSiC', '34:SPHEREY', '35:SPHEREJ', '36:SPHEREH', '37:SPHEREKs', '38:SPHEREJ2', '39:SPHEREJ3', '40:SPHEREH2', '41:SPHEREH3', '42:SPHEREK1', '43:SPHEREK2'))
    

    # interpolate linder models for 1 Neptune Mass planet to get limiting magnitude
    lines = np.where(np.logical_or(linder['2:Mass/Mearth'] == 10.0, linder['2:Mass/Mearth'] == 20.00))[0]
    lines = np.delete(lines, [28, 29, 30])

    filters = ['20:F200W', '22:F356W', '23:F444W', '24:F560W', '28:F1500W', '30:F2100W']
    filter_mags = drptst.interpolate_magnitude(linder_file, filters, lines)

    filter_mags = np.array(filter_mags)
    dist = np.array([5, 10, 20, 50, 100]) # distances in pc --> MIRI is better w/in 20pc (leave the further stuff to NIRCam)
    app_mags = np.array(drptst.absolute_to_apparent_magnitude(filter_mags, dist)) # the shape is (distances, ages, filters(magnitudes))

    """ 
    format of the different magnitudes in the grid: 
    -------------------------------------------------
    distances: [5, 10, 20, 50, 100] pc
    ages: [ 1,  1.2589,  1.5849,  1.9953,  2.5119, 3.1623,  3.9811,  5.0119,  6.3096,  7.9433, 10, 12.5893, 15.8489, 19.9526] in Myrs
    filters: ['20:F200W', '22:F356W', '23:F444W', '24:F560W', '28:F1500W', '30:F2100W']
    """

    if instrument == "NIRCam":
        # read in photometry
        mags_df = Table.read(BASE_DIR/'data_files'/'JADES'/'JADES_photometry.dat', format = 'ascii')

        # get rid of nan's and too faint galaxies
        too_faint = app_mags[:,:,2].max()
        mags_df['colour_circ3'] = mags_df['F200W_MAG_CIRC3'] - mags_df['F444W_MAG_CIRC3']
        colour_gals = np.array(mags_df['colour_circ3'])
        nan = np.where(np.logical_or(np.isnan(mags_df['colour_circ3']) == True, mags_df['F444W_MAG_CIRC3'] >= too_faint))[0]

        colours_gals = np.delete(colour_gals, nan)
        filt1_clipped = np.delete(mags_df['F200W_MAG_CIRC3'], nan) # nircam filt1 = f200w
        filt2_clipped = np.delete(mags_df['F444W_MAG_CIRC3'], nan) # nircam filt2 = f444w
        filt3_clipped = np.delete(mags_df['F356W_MAG_CIRC3'], nan) # nircam filt3 = f356w

        # number the galaxies for future reference
        mags_df['JADES ID'] = np.arange(0, len(mags_df['F200W_MAG_CIRC3'])) 
        gal_id_clipped = np.delete(mags_df['JADES ID'], nan)

        # get the coordinates for the galaxies
        ra = np.delete(np.array(mags_df['RA']), nan)
        dec = np.delete(np.array(mags_df['DEC']), nan)

        # get data in the three filters
        filt1_file = os.path.join(data_dir, "hlsp_jades_jwst_nircam_goods-s-deep_f200w_v2.0_drz.fits")
        filt2_file = os.path.join(data_dir, "hlsp_jades_jwst_nircam_goods-s-deep_f444w_v2.0_drz.fits")
        filt3_file = os.path.join(data_dir, "hlsp_jades_jwst_nircam_goods-s-deep_f356w_v2.0_drz.fits")

        data_filt1, wcs_filt1 = drptst.extract_data(filt1_file)
        data_filt2, wcs_filt2 = drptst.extract_data(filt2_file)
        data_filt3, wcs_filt3 = drptst.extract_data(filt3_file)

        # filter specific stuff
        filters = ['F356W', 'F200W', 'F444W']
        column_titles = "JADES_ID RA Dec FWHM_F356W FWHM_F200W FWHM_F444W F356W F200W F444W F200W-F444W"
        pix_to_as_filt23 = 0.063
        pix_to_as_filt1 = 0.031


    if instrument == "MIRI":
        # read in photometry 
        # mags_file = os.path.join(data_dir, "MIRI/SMILES_photometry.dat")
        mags_df = Table.read(BASE_DIR/'data_files'/'SMILES'/'SMILES_photometry.dat', format = 'ascii')

        # get rid of nan's and too faint galaxies
        too_faint = app_mags[:,:,4].max()
        mags_df['colour_kron'] = mags_df['F1500W_MAG_KRON'] - mags_df['F2100W_MAG_KRON']
        colour_gals = np.array(mags_df['colour_kron'])
        nan = np.where(np.logical_or(np.isnan(mags_df['colour_kron']) == True, mags_df['F1500W_MAG_KRON'] >= too_faint))[0]

        colours_gals = np.delete(colour_gals, nan)
        filt1_clipped = np.delete(mags_df['F1500W_MAG_KRON'], nan) # miri filt1 = f1500w
        filt2_clipped = np.delete(mags_df['F2100W_MAG_KRON'], nan) # miri filt2 = f2100w

        # number the galaxies for future reference
        mags_df['SMILES ID'] = np.arange(0, len(mags_df['F1500W_MAG_CIRC3']))
        gal_id_clipped = np.delete(mags_df['SMILES ID'], nan)

        # get the coordinates for the galaxies
        ra = np.delete(np.array(mags_df['RA']), nan)
        dec = np.delete(np.array(mags_df['DEC']), nan)

        # get data in the three filters
        filt1_file = os.path.join(data_dir, "hlsp_smiles_jwst_miri_goodss_f1500w_v1.0_drz.fits")
        filt2_file = os.path.join(data_dir, "hlsp_smiles_jwst_miri_goodss_f2100w_v1.0_drz.fits")
        filt3_file = os.path.join(data_dir, "hlsp_smiles_jwst_miri_goodss_f560w_v1.0_drz.fits")

        data_filt1, wcs_filt1 = drptst.extract_data(filt1_file)
        data_filt2, wcs_filt2 = drptst.extract_data(filt2_file)
        data_filt3, wcs_filt3 = drptst.extract_data(filt3_file)

        # filter specific stuff
        filters = ['F560W', 'F1500W', 'F2100W']
        column_titles = "SMILES_ID RA Dec FWHM_F560W FWHM_F1500W FWHM_F2100W F1500W F2100W F1500W-F2100W"
        pix_to_as_filt23 = 0.11
        pix_to_as_filt1 = 0.11


    n_gal = len(ra)

    # create empty arrays to store fwhms
    fwhms_filt3 = np.empty(n_gal)
    fwhms_filt1 = np.empty(n_gal)
    fwhms_filt2 = np.empty(n_gal)

    # empty arrays to store gaussian models
    gaussians_filt3 = np.empty(n_gal, dtype=object)
    gaussians_filt1 = np.empty(n_gal, dtype=object)
    gaussians_filt2 = np.empty(n_gal, dtype=object)

    dropout_check = np.empty(n_gal, dtype=object)

    # loop through all galaxies to check if they drop out
    for i, (ra_i, dec_i) in enumerate(zip(ra, dec)):

        fwhm_filt3, model_filt3 = drptst.get_fwhm(data_filt3, ra_i, dec_i, wcs_filt3)
        fwhm_filt1, model_filt1 = drptst.get_fwhm(data_filt1, ra_i, dec_i, wcs_filt1)
        fwhm_filt2, model_filt2 = drptst.get_fwhm(data_filt2, ra_i, dec_i, wcs_filt2)

        fwhms_filt3[i] = fwhm_filt3
        fwhms_filt1[i] = fwhm_filt1
        fwhms_filt2[i] = fwhm_filt2

        gaussians_filt3[i] = model_filt3
        gaussians_filt1[i] = model_filt1
        gaussians_filt2[i] = model_filt2

        dropout_check[i] = drptst.dropout_test([fwhm_filt3, fwhm_filt1, fwhm_filt2], [model_filt3, model_filt1, model_filt2], 2.5, 10)

    # these galaxies neeed to be checked by eye
    check_by_eye = np.where(dropout_check == "Check by eye")[0]
    print("Galaxies that need to be checked by eye:", check_by_eye)

    # plot each og these galaxy side-by-side in the three filters
    for i in check_by_eye:
        print(f"\n Galaxy {i}")
        drptst.plot_side_by_side(i, np.array([data_filt3, data_filt1, data_filt2]), ra[i], dec[i], np.array([wcs_filt3, wcs_filt1, wcs_filt2]), filters)

    # input which galaxies you want to KEEP
    user_input = input("\n Enter indicies of galaxies to KEEP (separated by commas): ")

    if user_input.strip() != "":
        keep_indices = np.array([int(x.strip()) for x in user_input.split(",")])

        for i in keep_indices:
            if i in check_by_eye:
                dropout_check[i] = "None"

    # keep galaxies that pass dropout check
    keep_sources = np.array(np.where(dropout_check == 'None')[0])

    # convert fwhm to arcsec
    fwhm_as_filt3 = fwhms_filt3*pix_to_as_filt23
    fwhm_as_filt1 = fwhms_filt1*pix_to_as_filt1
    fwhm_as_filt2 = fwhms_filt2*pix_to_as_filt23

    # save final photometry file for either instrument
    if instrument == "NIRCam":
        survey = 'JADES'
        save_photometry = np.column_stack((gal_id_clipped[keep_sources], ra[keep_sources], dec[keep_sources], fwhm_as_filt3[keep_sources], fwhm_as_filt1[keep_sources], 
                                fwhm_as_filt2[keep_sources], filt1_clipped[keep_sources], filt3_clipped[keep_sources], filt2_clipped[keep_sources], colours_gals[keep_sources]))

    if instrument == "MIRI":
        survey = 'SMILES'
        save_photometry = np.column_stack((gal_id_clipped[keep_sources], ra[keep_sources], dec[keep_sources], fwhm_as_filt3[keep_sources], fwhm_as_filt1[keep_sources], 
                                fwhm_as_filt2[keep_sources], filt1_clipped[keep_sources], filt2_clipped[keep_sources], colours_gals[keep_sources]))

    np.savetxt(BASE_DIR/'data_files'/f'{survey}'/f'{outputfilename}', save_photometry, header=column_titles,comments="", fmt="%s")
    print(f"Output written to '{outputfilename}'")

if __name__ == "__main__":
    main()
