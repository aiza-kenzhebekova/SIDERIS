""" 
Background object probability calculator functions
MPhys Project 
Aiza Kenzhebekova 
02/03/2026
"""
from astropy.table import Table
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import requests
import re
import copy
import time
from pathlib import Path


def run_trilegal_model(gc_l, gc_b, filters, field = '1', posturl = 'https://stev.oapd.inaf.it/cgi-bin/trilegal', resulturl = 'https://stev.oapd.inaf.it'): 
    """
    Run trilegal at a specific (l, b) and field of view. 
    Returns the URL for the data -- need to wait a bit for the link to be ready
    """
    # regex to get the result file name
    TRILEGAL_REGEX = re.compile(r'a href\S*.dat')

    TRILEGAL_DEFAULT_PARAMS = {
    'binary_frac': '0.3',
    'binary_kind': '1',
    'binary_mrinf': '0.7',
    'binary_mrsup': '1',
    'bulge_a': '1',
    'bulge_a0': '95',
    'bulge_am': '2500',
    'bulge_b': '-2.0e9',
    'bulge_csi': '0.31',
    'bulge_cutoffmass': '0.01',
    'bulge_eta': '0.68',
    'bulge_file': 'tab_sfr/file_sfr_bulge_zoccali_p03.dat',
    'bulge_kind': '2',
    'bulge_phi0': '15',
    'bulge_rho_central': '406.0',
    'eq_alpha': '0',
    'eq_delta': '0',
    'extinction_h_r': '100000',
    'extinction_h_z': '110',
    'extinction_infty': '0.0398',
    'extinction_kind': '2',
    'extinction_rho_sun': '0.00015',
    'extinction_sigma': '0',
    'field': '1',
    'gal_coord': '1',  # 1 -> set to gl_deg, gb_deg, 2 -> set to ra_hr, decl_deg
    'gc_b': '258.0',
    'gc_l': '-30.',
    'halo_a': '1',
    'halo_b': '0',
    'halo_file': 'tab_sfr/file_sfr_halo.dat',
    'halo_kind': '2',
    'halo_q': '0.65',
    'halo_r_eff': '2800',
    'halo_rho_sun': '0.00015',
    'icm_lim': '4',
    'imf_file': 'tab_imf/imf_chabrier_lognormal.dat',
    'mag_lim': '26',
    'mag_res': '0.1',
    'object_a': '1',
    'object_av': '1.504',
    'object_avkind': '1',
    'object_b': '0',
    'object_cutoffmass': '0.8',
    'object_dist': '1658',
    'object_file': 'tab_sfr/file_sfr_m4.dat',
    'object_kind': '0',
    'object_mass': '1280',
    'output_kind': '1',
    'photsys_file': 'tab_mag_odfnew/tab_mag_jwst_nircam_widemedium_nov22.dat',
    'r_sun': '8700',
    'submit_form': 'Submit',
    'thickdisk_a': '1',
    'thickdisk_b': '0',
    'thickdisk_file': 'tab_sfr/file_sfr_thickdisk.dat',
    'thickdisk_h_r': '2800',
    'thickdisk_h_z': '800',
    'thickdisk_kind': '0',
    'thickdisk_r_max': '15000',
    'thickdisk_r_min': '0',
    'thickdisk_rho_sun': '0.0015',
    'thindisk_a': '0.8',
    'thindisk_b': '0',
    'thindisk_file': 'tab_sfr/file_sfr_thindisk_mod.dat',
    'thindisk_h_r': '2800',
    'thindisk_h_z0': '95',
    'thindisk_hz_alpha': '1.6666',
    'thindisk_hz_tau0': '4400000000',
    'thindisk_kind': '3',
    'thindisk_r_max': '15000',
    'thindisk_r_min': '0',
    'thindisk_rho_sun': '59',
    'trilegal_version': '1.6',
    'z_sun': '24.2'
    }
   
    inputparams = copy.deepcopy(TRILEGAL_DEFAULT_PARAMS)

    inputparams['gc_l'] = str(gc_l)
    inputparams['gc_b'] = str(gc_b)
    inputparams['photsys_file'] = str(filters)
    inputparams['field'] = field

    req = requests.post(posturl, data=inputparams)
    resultfile = TRILEGAL_REGEX.search(req.text)

    if (resultfile):
        url = resulturl+resultfile[0].split('..')[1]
    else:
        print(req.text)
    
    return url


def download_trilegal_model(url, outfilename = 'trilegal_test.dat'):
    """
    Takes the URL from running trilegal and writes/saves the data file
    """
    response = requests.get(url, stream=True)
    
    with open(outfilename, mode="wb") as file:
        for chunk in response.iter_content(chunk_size=10 * 1024):
            file.write(chunk)

    return


def get_expectation_value(instrument, magnitude, colour, separation, trilegal_fov=2592000, trilegal_file="Default"):
    """ 
    Obtains the expectation value of an object with a certain colour, magnitude, and at a certain separation from the star 
    Works for NIRCam (F444W magnitude and F200W-F444W colour) and MIRI (F1500W magnitude and F1500W-F2100W colour)
    -----------------------------------------------------------------------------------------------------------------------
    Inputs:
    instrument = "MIRI" or "NIRCam"
    magnitude of the object 
    colour of the object
    separation of the object form the star (in arcsec)
    trilegal field of view (in sq. arcsec)
    """
    trilegal_fov = float(trilegal_fov)
    BASE_DIR = Path(__file__).resolve().parent.parent

    if instrument == 'NIRCam':
        # JADES data
        jades_mags = Table.read(BASE_DIR/'data_files'/'JADES'/'JADES_photometry_keep_default.dat', format = 'ascii') # change to your filename
        f200w = jades_mags['F200W']
        f444w = jades_mags['F444W']
        colour_jades = f200w - f444w
        exp_val_jades = len(np.where(np.logical_and(colour_jades > colour, f444w < magnitude))[0])

        # trilegal: NIRCam
        if trilegal_file=="Default":
            trilegal_nircam = Table.read(BASE_DIR/'data_files'/'trilegal_data'/'NIRCam_trilegal_53.11341_-27.92378_fov1296000.0.dat',
                     format='ascii') # change to your filename
        else:
            trilegal_nircam = Table.read(trilegal_file, format='ascii')

        f444w_trilegal = trilegal_nircam['F444W']
        f200w_trilegal = trilegal_nircam['F200W']
        colour_nircam_trilegal = f200w_trilegal - f444w_trilegal
        exp_val_trilegal_nircam = len(np.where(np.logical_and(colour_nircam_trilegal > colour, f444w_trilegal < magnitude))[0])
        # print(exp_val_trilegal_nircam)

        # JADES fov 68 armin^2 
        exp_val = (exp_val_jades * (np.pi * separation**2) / (68 * 3600)) + (exp_val_trilegal_nircam * (np.pi * separation**2) / (trilegal_fov))

    if instrument == 'MIRI':
        # SMILES data
        smiles_mags = Table.read(BASE_DIR/'data_files'/'SMILES'/'SMILES_photometry_keep_default.dat', format = 'ascii')  # change to your filename
        f1500w = smiles_mags['F1500W']
        f2100w = smiles_mags['F2100W']
        colour_smiles = smiles_mags['F1500W-F2100W']
        exp_val_smiles = len(np.where(np.logical_and(colour_smiles > colour, f1500w < magnitude))[0])

        # trilegal: MIRI
        if trilegal_file=="Default":
            trilegal_miri = Table.read(BASE_DIR/'data_files'/'trilegal_data'/'MIRI_trilegal_53.11644_-27.87666_fov1296000.0.dat',
                     format='ascii') # change to your filename
        else:    
            trilegal_miri = Table.read(trilegal_file, format='ascii')

        f1500w_trilegal = trilegal_miri['F1500W']
        f2100w_trilegal = trilegal_miri['F2100W']
        colour_miri_trilegal = f1500w_trilegal - f2100w_trilegal
        exp_val_trilegal_miri = len(np.where(np.logical_and(colour_miri_trilegal > colour, f1500w_trilegal < magnitude))[0])

        # SMILES fov is 34 arcmin^2
        exp_val = (exp_val_smiles * (np.pi * separation**2) / (34 * 3600)) + (exp_val_trilegal_miri * (np.pi * separation**2) / (trilegal_fov)) 

    return exp_val


def calculate_probability(input_file, trilegal_los=[3.3, 50], trilegal_fov=2592000, run_trilegal=False):
    """ 
    Calculates the probability of finding at least one object with specified colour, mag, 
        separation for a given line_of_sight (for trilegal stars)
    -------------------------------------------------------------------------------------
    Inputs: 
    input_file format is a .dat file where the header is either "F444W F200W-F444W Separation" or "F1500W F1500W-F2100W Separation" 
            depending on the instrument (separation in arcsec). This can run for multiple entries. 
        
    trilegal_los = line of sight for which you want to run trilegal (i.e RA, Dec of the target; [RA, Dec] format) 
    trilegal_fov = field of view for which you want to run trilegal (in sq. arcsec)

    Outputs: 
    Probabilities (percentage)
    """
    # read input file
    df = Table.read(input_file, format = 'ascii')

    # check if it's NIRCam or MIRI
    if df[0][:].name == 'F444W':
        instrument = 'NIRCam'
        magnitude = df['F444W']
        colour = df['F200W-F444W']
        separation = df['Separation']
        filters = 'tab_mag_odfnew/tab_mag_jwst_nircam_widemedium_nov22.dat'

    if df[0][:].name == 'F1500W':
        instrument = 'MIRI'
        magnitude = df['F1500W']
        colour = df['F1500W-F2100W']
        separation = df['Separation']
        filters = 'tab_mag_odfnew/tab_mag_jwst_miri_wide.dat'

    # run trilegal if needed
    if run_trilegal == True:
        # coords and fov for trilegal
        trilegal_fov_deg = trilegal_fov / 3600**2
        c_icrs = SkyCoord(ra=trilegal_los[0]*u.degree, dec=trilegal_los[1]*u.degree, frame='icrs')
        l = c_icrs.galactic.l.deg
        b = c_icrs.galactic.b.deg

        # run trilegal
        url = run_trilegal_model(l, b, filters, field=trilegal_fov_deg)
        print("Wait 2.5 mins for TRILEGAL")
        print(url)
        time.sleep(180)

        # save the data
        outfilename = f'data_files/trilegal_dat/{instrument}_trilegal_{trilegal_los[0]}_{trilegal_los[1]}_fov{trilegal_fov}.dat'
        download_trilegal_model(url, outfilename=outfilename)

    # w/o trilegal
    if run_trilegal == False:
        outfilename = "Default"

    probabilities = []

    for i in range(len(df)):
        # get probability for each row in the input file
        mag = magnitude[i]
        col= colour[i]
        sep = separation[i] 

        exp_val = get_expectation_value(instrument, mag, col, sep, trilegal_fov=trilegal_fov, trilegal_file=outfilename)
        prob = (1 - np.exp(-exp_val)) * 100
        probabilities.append(prob)

    return probabilities

