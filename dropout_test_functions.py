""" 
Drop Out Test Functions
MPhys Project 
Aiza Kenzhebekova 
02/03/2026
"""

from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from scipy.optimize import curve_fit
from astropy.visualization import simple_norm
from astropy.modeling import models, fitting

""" 
First two functions are to plot the Linder+2019 models and interpolate between the limiting magnitudes
"""

def interpolate_magnitude(model, filters, lines):
    """ 
    Interpolates the Linder+2019 models

    Filters are in this format: ['22:F356W', '20:F200W', '23:F444W', '24:F560W', '28:F1500W', '30:F2100W']
    """
 
    # read in linder model
    model = Table.read(model, format = 'ascii', names = ('1:log(Age/yr)', '2:Mass/Mearth', '3:Radius/Rjupiter', '4:Luminosity/Ljupiter', '5:Teff/K', '6:logg/cgs', '7:NACOJ', 
                '8:NACOH', '9:NACOKs', '10:NACOLp', '11:NACOMp', '12:CousinsR', '13:CousinsI', '14:WISE1', '15:WISE2', '16:WISE3', '17:WISE4', '18:F115W', '19:F150W', '20:F200W', 
                '21:F277W', '22:F356W', '23:F444W', '24:F560W', '25:F770W', '26:F1000W', '27:F1280W', '28:F1500W', '29:F1800W', '30:F2100W', '31:F2550W', '32:VISIRB87',
                '33:VISIRSiC', '34:SPHEREY', '35:SPHEREJ', '36:SPHEREH', '37:SPHEREKs', '38:SPHEREJ2', '39:SPHEREJ3', '40:SPHEREH2', '41:SPHEREH3', '42:SPHEREK1', '43:SPHEREK2'))
    
    model_new = model[lines]
    filter_mags = []
    for i in range(len(lines)): # interpolate
        interp_filter_mags = []
        if i % 2 == 0:
            for j in range(len(filters)):
                interp_filter_mags.append(((model_new[filters[j]][i] - model_new[filters[j]][i+1])*0.3) + model_new[filters[j]][i+1])
            filter_mags.append(interp_filter_mags)
        else:
            pass

    return filter_mags


def absolute_to_apparent_magnitude(abs_mags, dist):
    """ 
    Convert absolute to apparent magnitude
    """

    dmod = dmod = 5.*np.log10(dist / 10)
    apparent_mags = []
    for i in range(len(dist)):
        app_mags = abs_mags + dmod[i]
        apparent_mags.append(app_mags)

    return apparent_mags


def extract_data(filename):
    """ 
    This get the relevant data for both NIRCam and MIRI form the files.
    ---------------------------------------------------------------------
    Inputs: 
    Filename (fits file)

    Outputs: 
    data, wcs
    """

    hdu = fits.open(filename, format='ascii')
    data = hdu[1].data
    header = hdu[1].header
    wcs = WCS(header)

    return data, wcs


def get_cutout_data(data, ra, dec, wcs, size=(50,50)):
    """ 
    Get cutout data of any given size (in pixels; default 50x50) for the JADES/SMILES surveys
    ------------------------------------------------------------------------------------------
    data, wcs from extract_data
    ra, dec = cordinates of the galaxies form photometry files
    """
    sky_coord = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    cutout = Cutout2D(data, position=sky_coord, size=size, wcs=wcs)

    return cutout.data
    

def fit_gaussian(cutout_data):
    """ 
    fit a 2d gaussian to the data (already cut out)
    ------------------------------------------------
    returns the gaussian fit model (object), including: 
    amplitude, x_mean, y_mean, x_stddev, y_stddev, theta
    """

    data = np.array(cutout_data)
    data[np.isnan(data)] = 0 # set nan's to zero

    # create coordinate grid
    full_data = data
    height, width = full_data.shape
    x, y = np.meshgrid(np.arange(width), np.arange(height))
    cy, cx = height / 2, width / 2

    r = np.sqrt((x - cx)**2 + (y - cy)**2)

    # tuning parameter (brightness penalty per pixel of distance) --> helps it find sources in the centre
    alpha = 0.01 * np.nanstd(data)
    score = data - alpha * r

    # initial parameter guesses
    peak_guess = np.nanmax(full_data)
    y0_guess, x0_guess = np.unravel_index(np.nanargmax(score), full_data.shape)
    sigma_guess = 1

    g_init = models.Gaussian2D(
        amplitude=peak_guess,
        x_mean=x0_guess,
        y_mean=y0_guess,
        x_stddev=sigma_guess,
        y_stddev=sigma_guess,
        theta=0.0)
    
    g_init.x_mean.min = 0
    g_init.x_mean.max = width - 1 # -1 to avoid indeing errors 
    g_init.y_mean.min = 0
    g_init.y_mean.max = height - 1

    g_init.x_stddev.min = 0.5
    g_init.y_stddev.min = 0.5
    g_init.x_stddev.max = 0.5 * min(width, height)
    g_init.y_stddev.max = 0.5 * min(width, height)

    # fit model
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(
        g_init,
        x,
        y,
        full_data)

    return fitted_model


def get_fwhm(data, ra, dec, wcs, size=(50,50), plot_results=False):
    """ 
    Estimate fwhm of the cutout using fit_gaussian function
    ---------------------------------------------------------
    Inputs:
    data, wcs from extract_data
    ra, dec = cordinates of the galaxies form photometry files

    Outputs:
    fwhm (scalar)
    gaussian model (object)
    """
    # get data and fit a gaussian
    cutout_data = get_cutout_data(data, ra, dec, wcs)
    gaussian_model = fit_gaussian(cutout_data)

    # standard deviation of gaussian
    sigma_x = gaussian_model.x_stddev.value
    sigma_y = gaussian_model.y_stddev.value

    # geometric mean to preserve total area 
    sigma_eff = np.sqrt(sigma_x * sigma_y)
    
    # fwhm-gaussian standard deviation formula
    fwhm = 2.355 * sigma_eff   

    # can plot results if desired
    if plot_results == True:

        ny, nx = cutout_data.shape
        y, x = np.mgrid[:ny, :nx]
        x0 = gaussian_model.x_mean.value
        y0 = gaussian_model.y_mean.value

        radius = 0.5 * fwhm
        theta = np.linspace(0, 2*np.pi, 200)

        x_circle = x0 + radius * np.cos(theta)
        y_circle = y0 + radius * np.sin(theta)

        model_image = gaussian_model(x, y)
        residual = cutout_data - model_image

        # colourbar
        vmin = np.nanmin(cutout_data)
        vmax = np.nanmax(cutout_data)
        norm = simple_norm(cutout_data, 'sqrt', min_cut=vmin, max_cut=vmax)

        # plot original cutout, gaussian fit, and residual
        fig, axes = plt.subplots(1, 3, figsize=(10, 4))
        axes[0].imshow(cutout_data, origin="lower", norm=norm, cmap='magma')
        axes[0].plot(x_circle, y_circle, color="red", linewidth=2, label="FWHM$_{eff}$/2")
        axes[1].imshow(model_image, origin="lower", norm=norm, cmap='magma')
        axes[1].plot(x_circle, y_circle, color="red", linewidth=2, label="FWHM$_{eff}$/2")
        axes[2].imshow(residual, origin="lower", norm=norm, cmap='magma')
        axes[2].plot(x_circle, y_circle, color="red", linewidth=2, label="FWHM$_{eff}$/2")

    return fwhm, gaussian_model


def dropout_test(fwhm_3, gaussian_models_3, max_fwhm_ratio, max_offset, size=(50,50)):
    """ 
    Tests if source drops out
    --------------------------
    Inputs: 
    fwhm (set of tree, one per filter)
    gaussian models (set of threee, one per filter)
    maximum fwhm ratio: the ratio between the forced photometry filter fwhm and other filter at which 
                            a source is deemed to drop out
    maximum offset difference: max distance between the source in the forced photometry filter and other 
                            filter before a source is deemed to drop out

    Outputs:
    Whether and how a source drops out. Can also return "Check by eye" for fringe cases.
    """
    # fwhm ratio
    fwhm_ratio_filt1 = fwhm_3[1]/fwhm_3[0] 
    fwhm_ratio_filt2 = fwhm_3[2]/fwhm_3[0]
    fwhm_condition = np.logical_or(fwhm_ratio_filt1 > max_fwhm_ratio, fwhm_ratio_filt2 > max_fwhm_ratio)

    # offset
    x0 = gaussian_models_3[0].x_mean.value
    y0 = gaussian_models_3[0].y_mean.value

    x1 = gaussian_models_3[1].x_mean.value
    y1 = gaussian_models_3[1].y_mean.value

    x2 = gaussian_models_3[2].x_mean.value
    y2 = gaussian_models_3[2].y_mean.value

    dx1 = x1 - x0
    dy1 = y1 - y0
    offset1 = np.hypot(dx1, dy1)

    dx2 = x2 - x0
    dy2 = y2 - y0
    offset2 = np.hypot(dx2, dy2)

    offset_condition = np.logical_or(offset1 > max_offset, offset2 > max_offset)

    # check by eye
    ny, nx = size
    y, x = np.mgrid[:ny, :nx]

    model_image_filt3 = gaussian_models_3[0](x, y)
    model_image_filt1 = gaussian_models_3[1](x, y)
    model_image_filt2 = gaussian_models_3[2](x, y)

    difference_filt3 = np.nanmax(model_image_filt3) - np.nanmin(model_image_filt3)
    difference_filt1 = np.nanmax(model_image_filt1) - np.nanmin(model_image_filt1)
    difference_filt2 = np.nanmax(model_image_filt2) - np.nanmin(model_image_filt2)

    # flaged if source is faint in drop out filter and not in other filters
    gaussian_sanity_check = np.logical_or((1.5 * difference_filt3) < difference_filt1, (1.5 * difference_filt3) < difference_filt2)

    # dropout evaluations
    if np.logical_and(offset_condition == True, fwhm_condition == True):
        return 'Drop Out--Offset and FWHM'
    
    elif offset_condition == True:
        return 'Drop Out--Offset only'
    
    elif fwhm_condition == True:
        if gaussian_sanity_check == True:
            return 'Check by eye'
        else:
            return 'Drop Out--FWHM only'
    
    else:
        return 'None'


def plot_side_by_side(i, data_3, ra, dec, wcs_3, filt, window_size=50, size=(50,50)):
    """ 
    Plots a source side-by-side in all three filters
    -------------------------------------------------
    Inputs: 
    i = index of the source
    data_3, wcs_3 = data and wcs's in all three filters from extract data
    ra, dec from photometry

    Outputs: 
    1x3 plot
    """
    # get cutout data in all three filters
    cutout_data_filt3 = get_cutout_data(data_3[0], ra, dec, wcs_3[0], size=size)
    fwhm_filt3, gaussian_model_filt3 = get_fwhm(data_3[0], ra, dec, wcs_3[0], size=size)

    cutout_data_filt1 = get_cutout_data(data_3[1], ra, dec, wcs_3[1], size=size)
    fwhm_filt1, gaussian_model_filt1 = get_fwhm(data_3[1], ra, dec, wcs_3[1], size=size)

    cutout_data_filt2 = get_cutout_data(data_3[2], ra, dec, wcs_3[2], size=size)
    fwhm_filt2, gaussian_model_filt2= get_fwhm(data_3[2], ra, dec, wcs_3[2], size=size)

    # find centre
    x0_filt3 = gaussian_model_filt3.x_mean.value
    y0_filt3 = gaussian_model_filt3.y_mean.value

    x0_filt1 = gaussian_model_filt1.x_mean.value
    y0_filt1 = gaussian_model_filt1.y_mean.value

    x0_filt2 = gaussian_model_filt2.x_mean.value
    y0_filt2 = gaussian_model_filt2.y_mean.value   

    # draw circle that corresponds to the fwhm 
    theta = np.linspace(0, 2*np.pi, 200)
    radius_filt3 = 0.5 * fwhm_filt3
    x_circle_filt3 = x0_filt3 + radius_filt3 * np.cos(theta)
    y_circle_filt3 = y0_filt3 + radius_filt3 * np.sin(theta)

    radius_filt1 = 0.5 * fwhm_filt1
    x_circle_filt1 = x0_filt1 + radius_filt1 * np.cos(theta)
    y_circle_filt1 = y0_filt1 + radius_filt1 * np.sin(theta)

    radius_filt2 = 0.5 * fwhm_filt2
    x_circle_filt2 = x0_filt2 + radius_filt2 * np.cos(theta)
    y_circle_filt2 = y0_filt2 + radius_filt2 * np.sin(theta)

    # colourbar
    norm_data_filt3 = simple_norm(cutout_data_filt3, 'sqrt', min_cut=np.nanmin(cutout_data_filt3), max_cut=np.nanmax(cutout_data_filt3))
    norm_data_filt1 = simple_norm(cutout_data_filt1, 'sqrt', min_cut=np.nanmin(cutout_data_filt1), max_cut=np.nanmax(cutout_data_filt1))
    norm_data_filt2 = simple_norm(cutout_data_filt2, 'sqrt', min_cut=np.nanmin(cutout_data_filt2), max_cut=np.nanmax(cutout_data_filt2))
    
    # plot
    fig, axes = plt.subplots(1, 3, figsize=(8, 5))

    axes[0].imshow(cutout_data_filt3, origin="lower", norm=norm_data_filt3, cmap='magma')
    axes[0].plot(x_circle_filt3, y_circle_filt3, color="red", linewidth=1, label="FWHM")
    axes[0].axvline(x0_filt3, color='black', linestyle='dashed', linewidth=0.5)
    axes[0].axhline(y0_filt3, color='black', linestyle='dashed', linewidth=0.5)
    axes[0].set_title(f'{i}: {filt[0]}')

    axes[1].imshow(cutout_data_filt1, origin="lower", norm=norm_data_filt1, cmap='magma')
    axes[1].plot(x_circle_filt1, y_circle_filt1, color="red", linewidth=1, label="FWHM")
    axes[1].axvline(x0_filt1, color='black', linestyle='dashed', linewidth=0.5)
    axes[1].axhline(y0_filt1, color='black', linestyle='dashed', linewidth=0.5)
    axes[1].set_title(f'{filt[1]}')

    axes[2].imshow(cutout_data_filt2, origin="lower", norm=norm_data_filt2, cmap='magma')
    axes[2].plot(x_circle_filt2, y_circle_filt2, color="red", linewidth=1, label="FWHM")
    axes[2].axvline(x0_filt2, color='black', linestyle='dashed', linewidth=0.5)
    axes[2].axhline(y0_filt2, color='black', linestyle='dashed', linewidth=0.5)
    axes[2].set_title(f'{filt[2]}')

    plt.show()