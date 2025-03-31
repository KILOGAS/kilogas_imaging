"""
Imaging pipeline for KILOGAS data.
Tested on CASA 6.5.4 and MSs from ALMA Cycle 11.
"""

import configparser
import glob
import json
import os
import re
import shutil
import sys

import analysisUtils as au
import astropy.constants as const
import astropy.units as u
import numpy as np

from astropy.io import fits
from astropy.table import Table, Column
from datetime import datetime

# ==============================================================================

def parse_config(config):
    """
    Reads an input configuration file name (e.g. 'config.ini'), extracts
    parameter values and returns list of parameters values for each section of
    the config file.

    Sections:
    - [paths]: Main directory and subdirectory managemaent
    - [input]: Select targets and bands to image
    - [imaging]: A subset of parameters and options for tclean()
    - [output]: Set various output options
    """
    # ---=== Initialise ConfigParser ===---
    cfg = configparser.ConfigParser(inline_comment_prefixes='#')
    # Read configuration file
    cfg.read(config)

    # ---=== Get items from each section ===---
    # [paths]:
    # Path to top-level directory of ALMA data (usually the project code)
    d_pth = cfg.get('paths', 'data_path')
    # Path to top-level directory of imaging products
    p_pth = cfg.get('paths', 'prod_path')
    # Name to give top-level directory of output (e.g. DR1, test_01)
    o_dir = cfg.get('paths', 'dir_name')
    if o_dir == 'auto':
        # Create unique name based on current date and time
        pipeline_start = datetime.utcnow().strftime('%Y-%m-%d-%H%M%S')
        o_dir = 'kgas_pipeline_' + pipeline_start
        '''
        NOTE: This option is good for testing short runs of the pipeline with
              subtle changes in parameters, as it will prevent overwriting the
              directory.
        '''
    else:
        # Remove '/' from name in case one was added by user
        o_dir = o_dir.strip('/')

    # [input]:
    # Target selection string
    trgts = cfg.get('input', 'targets')
    if trgts=='all':
        # Image all available targets with 12m MSs in table
        trgts = 'all'
    elif trgts.endswith('.txt'):
        # Read .txt file containing list of targets, convert to list
        sources_file = np.loadtxt(trgts, dtype=int, delimiter=',')
        trgts = [int(id) for id in sources_file]
        '''
        NOTE: This is a comma delimited list of integers in a plain text file.
        '''
    else:
        # Convert string to list of targets ids (integers)
        trgts = [int(id) for id in re.split(', |,| ', trgts.strip('()[]{}'))]
    # Get list of which bands to image as cubes
    lbstr = cfg.get('input', 'line_bands').strip('()[]{}')
    if lbstr == 'none':
        l_bnd = [] # If the user only wants to image continuum
    else:
        l_bnd = re.split(', |,| ', lbstr)
    # Get list of which bands to jointly image as continuum using MFS
    c_bnd = re.split(', |,| ', cfg.get('input', 'cont_bands').strip('()[]{}'))
    '''
    NOTE: These line/band short names need to be consistent with the dictionary
          of transition/band labels.
    '''

    # [imaging]:
    # Read the cell size in arcseconds
    cl_sz = cfg.getfloat('imaging', 'cell_size')
    # Noise threshold for stopping tclean(), in multiples of the measured RMS
    n_rms = cfg.getfloat('imaging', 'rms_threshold')
    # Get channel width for imaging, in km/s
    chn_w = cfg.getfloat('imaging', 'chan_width')
    # Number of channels to use for RMS
    l_nch = re.split(', |,| ', cfg.get('imaging', 'nchan_rms').strip('()[]{}'))
    # Format into tuple, corresponding to number of channels from each end
    l_nch = int(l_nch[0]), int(l_nch[-1])
    # Weighting scheme for line imaging (e.g. briggs)
    ln_wt = cfg.get('imaging', 'line_weight')
    # Weighting scheme for continuum imaging (e.g. natural)
    cn_wt = cfg.get('imaging', 'cont_weight')
    # Robust parameter for when briggs weighting is used
    b_rob = cfg.getfloat('imaging', 'robust')
    # Limit max scale used by asp. This will be set to the MRS
    lim_s = cfg.getboolean('imaging', 'limit_scales')
    # Allow the pipeline to perform imcontsub() where necessary
    c_sub = cfg.getboolean('imaging', 'do_contsub')

    # [output]:
    # Export FITS files of imaging products?
    x_fts = cfg.getboolean('output', 'export_fits')
    # List of image types to export as FITS. If default then use predefined list
    l_fts = re.split(', |,| ', cfg.get('output', 'fits_list').strip('()[]{}'))
    if l_fts == ['default']:
        l_fts = ['.image', '.image.pbcor', '.residual',
                 '.pb', '.model', '.mask']
    else:
        l_fts = l_fts
    # Convert cubes from Jy/beam -> K?
    cnv_k = cfg.getboolean('output', 'cube_in_K')
    # Use optical velocity convention in FITS files?
    opt_f = cfg.getboolean('output', 'use_optical')
    # Bunch up variables into list for easier handling later
    cfg_paths = [d_pth, p_pth, o_dir]
    cfg_input = [trgts, l_bnd, c_bnd]
    cfg_imaging = [cl_sz, n_rms, chn_w, l_nch,
                   ln_wt, cn_wt, b_rob, lim_s, c_sub]
    cfg_output = [x_fts, l_fts, cnv_k, opt_f]
    '''
    NOTE: This config parser is a bit rough, but it will do for now. Some
          things to improve would be:
              1. Add a way to handle incorrect/invalid options before they
                 reach the point where they are used in the rest of the code.
              2. Use dictionaries as opposed to lists so that extracting
                 options isn't reliant on list indexing (not best practice, but
                 this part of the pipeline ended up being a bit rushed).
              3. Handle more of the variable setting and interpreting in this
                 function, so that it returns more ready-to-use variables.
    '''

    return cfg_paths, cfg_input, cfg_imaging, cfg_output


def unpack_path(pathstr, delimiter='|'):
    """
    Splits concatenated string of multiple paths, delimited by '|'.
    """
    path_list = pathstr.split(delimiter)
    # Return empty list if no path found
    if path_list == ['']:
        path_list = []

    return path_list


def MRS_L5(vis, field, spw):
    """
    Returns the maximum resolvable scale (in arcseconds) using the 5th
    percentile of u,v distances (L_5) in metres. This is an empirically derived
    formula from the ALMA Cycle 11 Technical Handbook (eq. 7.6):

    theta_MRS = 0.983 * wavelength / L_5 [in radians]

    Input a path to a measurement set, the field name, and desired spw.
    - The field name corresponds to the target's name within the MS.
    - The spw argument is a string, and follows the same spw/channel
    format as tclean(), i.e. spw 25, channels [800, 1200): spw='25:800~1199'
    """
    # Get baseline statistics
    L_stats = au.getBaselineStats(msFile=vis, percentile=5, field=field)
    # Get 5th percentile of baseline lengths
    L_5 = L_stats[0] * u.m

    # Get mean frequency of spw
    if ',' in spw:
        mean_f = au.getMeanFreqOfSpwlist(vis=vis, spw=spw) * u.Hz
    else:
        mean_f = au.getMeanFreqOfSpwSelection(vis=vis, spw=spw) * u.Hz
    '''
    NOTE: Use au.getMeanFreqOfSpwlist() when more than one spw is given. It
          is possible to use au.getMeanFreqOfSpwSelection(), and it would be
          better to use it, but some reformatting of the string is required.
    '''

    # Convert frequency to wavelength
    wl = mean_f.to(u.m, equivalencies=u.spectral())

    # Calculate MRS (and convert to arcseconds)
    theta_MRS = (((0.983 * wl / L_5) * u.rad).to(u.arcsec)).value

    return theta_MRS


# ==============================================================================

# ---=== Initialise pipeline ===---
# Get location of repository
rpath = os.path.abspath(__file__).strip('kilogas_imaging.py')

# Parse config file and get parameters
cfg_lists = parse_config(config=rpath + 'config.ini')
# Get section options from parsed config file, and set appropriate variables
cfg_paths, cfg_input, cfg_imaging, cfg_output = cfg_lists

# Path/directory options
dpath = cfg_paths[0] # Absolute path to ALMA project folder
ppath = cfg_paths[1] # Absolute path to desired product location
dname = cfg_paths[2] # Name of directory to house all imaging products

# Target/line selection options
source_list = cfg_input[0] # List of sources to image, or choose all available
line_bands = cfg_input[1] # List of bands to image as cubes
cont_bands = cfg_input[2] # List of bands to jointly image as continuum

# Imaging options
c_size = cfg_imaging[0] # Cell size in arcseconds
rms_factor = cfg_imaging[1] # Mutliple of the measured RMS for clean threshold
ch_width = cfg_imaging[2] # Channel width in km/s
l_Nch, h_Nch = cfg_imaging[3] # Number of channels from left/right end for RMS
l_weighting = cfg_imaging[4] # Weighting scheme for line imaging
c_weighting = cfg_imaging[5] # Weighting scheme for continuum imaging
robust = cfg_imaging[6] # Value of robust parameter when using briggs weighting
limit_scales = cfg_imaging[7] # Flag for whether to limit max asp scale
do_contsub = cfg_imaging[8] # Flag for whether to allow running imconstub()

# FITS export options
export_fits = cfg_output[0] # Flag to export FITS of imaging products
export_list = cfg_output[1] # List of image types to export as FITS
cube_in_K = cfg_output[2] # Flag to do Jy/beam -> K conversion
use_optical = cfg_output[3] # Flag to use optical velocity in FITS files


# Import table of KILOGAS sources and paths to their MSs
kgas_tab = Table.read(rpath + 'tables/kilogas_ms_paths.fits', format='fits')

# Import dictionary of short names for each band/line
with open(rpath + 'band_names.json', 'r') as file:
    transitions_dict = json.loads(file.read())

# Import dictionary of automasking paramteres
with open(rpath + 'automask_params.json', 'r') as file:
    automask_params = json.loads(file.read())


# Max number of iterations for tclean()
n_iter = 10000000 # Large number so tclean() hits noise threshold not iterations

# ---=== Initialise summary tables ===---
# Base summary table to build other tables from
base_sum_tab = kgas_tab['KGAS_ID', 'IAUname', 'SB_num',
                        'has12m?', '12m_path', 'has7m?', '7m_path']
# Get length of summary table
bstl = len(base_sum_tab)
# Initialise columns
col_list = [Column(name='TA_str', length=bstl, dtype='S6',
                   description='Array setup as a string (12m/7m+12m)'),
            Column(name='beam_major', length=bstl, dtype=float, unit=u.arcsec,
                   description='Beam major axis'),
            Column(name='beam_minor', length=bstl, dtype=float, unit=u.arcsec,
                   description='Beam minor axis'),
            Column(name='beam_PA', length=bstl, dtype=float, unit=u.degree,
                   description='Beam position angle'),
            Column(name='RMS', length=bstl, dtype=float, unit=u.Jy / u.beam,
                   description='Measured root mean square noise'),
            Column(name='peak_SNR', length=bstl, dtype=float,
                   description='Peak signal to noise value'),
            Column(name='low_SNR_flag', length=bstl, dtype=bool,
                   description='True if peak_SNR < 5'),
            Column(name='MRS_L5', length=bstl, dtype=float, unit=u.arcsec,
                   description='Maximum resolvable scale'),
            Column(name='cont_bands', length=bstl, dtype='S128',
                   description='Bands used for continuum, "|" delimited'),
            Column(name='cont_sub_flag', length=bstl, dtype=bool,
                   description='True if continuum subtraction was performed'),
            Column(name='empty_model_flag', length=bstl, dtype=bool,
                   description='True if tclean() model is empty'),
            Column(name='baseline_min', length=bstl, dtype=float, unit=u.m,
                   description='Minimum projected baseline length'),
            Column(name='baseline_max', length=bstl, dtype=float, unit=u.m,
                   description='Maximum projected baseline length'),
            Column(name='baseline_5th', length=bstl, dtype=float, unit=u.m,
                   description='5th percentile projected baseline length'),
            Column(name='baseline_80th', length=bstl, dtype=float, unit=u.m,
                   description='80th percentile projected baseline length'),
            Column(name='flux_calibrator', length=bstl, dtype='S128',
                   description='Flux calibators, "|" delimited'),
            Column(name='phase_calibrator', length=bstl, dtype='S128',
                   description='Phase calibrators, "|" delimited'),
            Column(name='bandpass_calibrator', length=bstl, dtype='S128',
                   description='Bandpass calibrators, "|" delimited'),
            Column(name='imaged_flag', length=bstl, dtype=bool,
                   description='True if target has been imaged')]
# Add columns to base summary table
base_sum_tab.add_columns(col_list)
'''
NOTE: So this is far from elegant and I hate it, but the summary tables need
      to be initalised somehow.
'''
# Set names for summary tables
cont_tab_name = dname + '_continuum_table.fits' # Continuum summary table
line_tab_names = [] # Line summary tables
for l_band in line_bands:
    line_tab_names.append(dname + '_' + l_band + '_' +
                          str(ch_width) + 'kmps' + '_table.fits')

# Check if tables exist in output dir
if os.path.exists(ppath + dname + '/original/' + cont_tab_name):
    # Import continuum summary table
    c_sum_tab = Table.read(ppath + dname + '/original/' + cont_tab_name)
else:
    # Initialise continuum summary table
    c_sum_tab = base_sum_tab.copy()
    c_sum_tab.remove_column('cont_sub_flag') # Remove only for continuumm

l_sum_tab = [] # List of line summary tables
for tn in line_tab_names:
    if os.path.exists(ppath + dname + '/original/' + tn):
        # Import line summary table (for each band)
        l_sum_tab.append(Table.read(ppath + dname + '/original/' + tn))
    else:
        # Initialise line summary table (for each band)
        ls_tab = base_sum_tab.copy()
        ls_tab.remove_column('cont_bands') # Remove only for line tables
        l_sum_tab.append(ls_tab)
'''
NOTE: This is for when you may want to re-run the pipeline for a few sources,
      and hence only want to update the table and not overwrite it completely.
'''

# ---=== Select sources to image from main table ===---
# Check which targets have 12m MSs
all_kgas_ids = kgas_tab['KGAS_ID'].value # Get all KGAS IDs
all_12m_flag = kgas_tab['has12m?'].value # Flag for if 12m MS is found

if source_list == 'all':
    # Get list of all available targets with 12m MSs
    source_list = all_kgas_ids[all_12m_flag].tolist() # Update source_list
    print('Imaging all KGAS targets that have 12m MSs.')
    no_ms_list = all_kgas_ids[~all_12m_flag].tolist()
    if len(no_ms_list) > 0:
        print('No 12m MSs found for KGAS {}.\n'.format(no_ms_list) +
              'These targets will not be imaged.\n')
elif type(source_list) == list:
    # Check the table and remove any sources in this list without 12m MSs
    ms_avail = np.isin(source_list, all_kgas_ids[all_12m_flag])
    no_ms_list = np.array(source_list)[~ms_avail].tolist()
    source_list = np.array(source_list)[ms_avail].tolist() # Update source_list
    print('Imaging KGAS {}.'.format(source_list))
    if len(no_ms_list) > 0:
        print('No 12m MSs found for KGAS {}.\n'.format(no_ms_list) +
              'These targets will not be imaged.\n')

# Create subset of imported table using those listed in source_list
tab_mask = np.isin(kgas_tab['KGAS_ID'], source_list)
# Create subset of table for imaging using mask
im_tab = kgas_tab[tab_mask]
# Get columns from table
kgas_id = im_tab['KGAS_ID']
iau_name = im_tab['IAUname']
has_12m = im_tab['has12m?']
path_12m = Column(np.ma.getdata(im_tab['12m_path']))
has_7m = im_tab['has7m?']
path_7m = Column(np.ma.getdata(im_tab['7m_path']))

# ---=== Begin pipeline ===---
# Iterate over all sources to be imaged
for i, kg_id in enumerate(kgas_id):
    # ---=== Create string of KILOGAS ID ===---
    g_name = 'KGAS' + str(kg_id)

    # ---=== Create mask for continuum summary tables ===---
    stm = c_sum_tab['KGAS_ID'] == kg_id # Mask should also work for line tables

    # ---=== Check for array configuration ===---
    # Set t_arr for filenames etc
    if has_12m[i] & ~has_7m[i]:
        t_arr = '12m'
    elif has_12m[i] & has_7m[i]:
        t_arr = '7m+12m'
    # Add to summary tables
    c_sum_tab['TA_str'][stm] = t_arr
    for ti in range(len(l_sum_tab)):
        l_sum_tab[ti]['TA_str'][stm] = t_arr

    # ---=== Create directory structure, remove conflicting files ===---
    dir_name = dname + '/' + 'original' + '/' + g_name + '/'

    if os.path.exists(ppath + dir_name):
        print('Directory {} already exists.\n'.format(dir_name) +
              'Deleting any existing conflicting imaging products.\n')
        lbs = line_bands + ['continuum']
        for lb in lbs:
            # Get list of all possible conflicting files
            to_delete = glob.glob(ppath + dir_name + g_name + '_' + lb + '_*')
            # All previously imaged cubes
            prev_vw = glob.glob(ppath + dir_name + g_name + '_' + lb +
                                   '_' + '*kmps' + '_*')
            # Imaged cubes with the same velocity width
            curr_vw = glob.glob(ppath + dir_name + g_name + '_' + lb +
                                   '_' + str(ch_width) + 'kmps' + '_*')
            if len(curr_vw) > 0:
                # Remove previously imaged cubes with same channel width
                to_delete = list(set(to_delete) -
                                 (set(prev_vw) - set(curr_vw)))
                '''
                NOTE: This is very contrived. The aim here is to only delete
                      previously imaged cubes with the same velocity resolution
                      as the current run, so that you can run the pipeline
                      multiple times with different channel widths and compare
                      the outputs.
                '''
            if len(to_delete) > 0:
                for td in to_delete:
                    if os.path.isfile(td):
                        os.remove(td)
                    elif os.path.isdir(td):
                        shutil.rmtree(td)
        '''
        NOTE: Deleting previous all previous imaging products is an
              over-the-top safety measure to avoid potential later confusion.
              The only files that cause major issues by existing already are
              ".residual" and ".model", because tclean() will continue from
              these, effectively resuming the previous run of tclean().
        '''
    else:
        print('Creating directory {}\n'.format(dir_name))
        os.makedirs(ppath + dir_name)

    # ---=== Calculate average field position for phasecenter ===---
    # Get list of paths to all MSs for current target
    g_ms = unpack_path(path_12m[i]) + unpack_path(path_7m[i])
    # Flag for MS type (this is convoluted, I know)
    g_ms_flag = (['12m'] * len(unpack_path(path_12m[i])) +
                 ['7m'] * len(unpack_path(path_7m[i])))
    # Get list of paths to 12m MSs for target (useful for some functions)
    g_ms_12m = unpack_path(path_12m[i])
    # Iterate over each 12m MS path to get info from each MS for phasecentre
    dDir_list = []  # List of target field centre coords for each MS
    for mspath in g_ms_12m:
        # Get field indices for given source name from MS
        pc_tool = au.createCasaTool(msmdtool) # Initialises msmd tool
        pc_tool.open(mspath)
        src_ind = pc_tool.fieldsforname(iau_name[i])
        pc_tool.close() # Close msmd tool
        # Get RA and Dec from table
        tb_tool = au.createCasaTool(tbtool) # Initialises table tool
        tb_tool.open(mspath+ '/FIELD')
        delayDir = tb_tool.getcol('DELAY_DIR') # [RA/Dec, N fields, M sources]
        dDir_list.append(delayDir[:, :, src_ind]) # Add to list of coords
        tb_tool.close() # Close table tool
    # Concat coord arrays for phase centre calculation
    dDir_arr = np.concatenate(dDir_list, axis=2)
    # Calculate RA and Dec
    ra = dDir_arr[0, 0] * (12 / np.pi)  # RA in radians -> decimal HMS
    for j in range(len(ra)):
        if ra[j] < 0:
            ra[j] += 24 # Catch for if RA in decimal HMS is negative
    ra *= 15  # Convert RA in decimal HMS -> degrees
    dec = dDir_arr[1, 0] * (180 / np.pi)  # Dec in radians -> degrees
    # Calculate average RA and Dec for fields (only useful for mosaics)
    ra_avg = np.mean(ra)
    dec_avg = np.mean(dec)
    # Format calculated phase center to string
    phase_center = 'J2000 ' + str(ra_avg) + 'deg ' + str(dec_avg) + 'deg'

    # ---=== Calculate optimal image size ===---
    # Check 12m MSs for which spw is the lowest frequency
    min_12m_avf = [] # Minimum average frequency of all spws, per MS
    min_12m_spw = [] # Corresponding spw for each MS
    for mspath in g_ms_12m:
        t_dict = au.transitions(mspath, source=iau_name[i],
                                intent='OBSERVE_TARGET#ON_SOURCE')
        spws = list(t_dict) # List of spw IDs
        f_list = [] # List of mean frequencies
        for sw in spws:
            f_list.append(np.mean(au.getFrequencies(inputMs=mspath, spwId=sw)))
        # Get index of minimum frequency
        min_12m_avf.append(f_list[np.argmin(f_list)])
        min_12m_spw.append(spws[np.argmin(f_list)])
    # Get index of lowest freq. spw to use for optimal image size calculation
    lowest_f_ind = np.argmin(min_12m_avf)
    imsize_mspath = g_ms_12m[lowest_f_ind] # MS to use
    imsize_spw = min_12m_spw[lowest_f_ind] # Which spw to use

    # Calculate optimal image size
    cell_size = str(c_size)+'arcsec' # Set cell size, formatted for tclean()
    # Estimate extent of FoV in RA and Dec
    pltmos = au.plotmosaic(vis=imsize_mspath,
                           sourceid=iau_name[i],
                           spw=imsize_spw,
                           pblevel=0.05,
                           doplot=False,
                           verbose=False)
    '''
    NOTE: pblevel=0.05 allows for a good buffer around the imaged FoV
          (in tclean() the PB limit will be set to 0.2). Shrinking the image
          area would lead to warping of the FoV due to periodic wrapping.
    '''
    ra_max, ra_min, dec_max, dec_min = pltmos[1:]
    ra_FoV = ra_max - ra_min
    dec_FoV = dec_max - dec_min
    # Calculate optimal image size
    opt_ra = au.getOptimumSize(int(np.ceil(abs(ra_FoV)/c_size)))
    opt_dec = au.getOptimumSize(int(np.ceil(abs(dec_FoV)/c_size)))
    im_size = [opt_ra, opt_dec]
    '''
    NOTE: This is being done here so that all images for all bands are on
          the same grid/array for easier comparison. The spw with the lowest
          frequency is being chosen to maximise the image size, as the FoV
          is inversely proportional to frequency.

          If you want to be able to sample all beams with an equal number of
          pixels (i.e. variable cell size), I can add this functionality
          relatively easily (upon request).
    '''
    # ---=== Get baseline statistics for summary tables ===---
    # Get list of all baseline lengths for all MSs for the current target
    bl_list = au.getBaselineLengthsMultiVis(vislist=g_ms, field=iau_name[i])
    '''
    NOTE: One concern with this function however, is how the statistics are
          skewed when there are multiple MSs for each array. Some further
          testing is required, so don't take these baseline stats as being
          perfect. The projected baselines won't be perfect anyway I think, as
          this function calculates the projected baseline lengths for one scan.

          This function would have been very helpful in other parts of the
          code, but I found it too late in development. An improvement for the future would be to replace instances where au.getBaselineStats() is
          iterated over multiple times with different MSs, as it takes a long
          time to run each time.
    '''
    # Get minimum, 5th percentile, 80th percentile, maximum baseline lengths
    bl_pcnt = np.nanpercentile(bl_list, (0, 5, 80, 100)) * u.m
    # Add stats to summary tables
    c_sum_tab['baseline_min'][stm] = bl_pcnt[0]
    c_sum_tab['baseline_max'][stm] = bl_pcnt[3]
    c_sum_tab['baseline_5th'][stm] = bl_pcnt[1]
    c_sum_tab['baseline_80th'][stm] = bl_pcnt[2]
    for ti in range(len(l_sum_tab)):
        l_sum_tab[ti]['baseline_min'][stm] = bl_pcnt[0]
        l_sum_tab[ti]['baseline_max'][stm] = bl_pcnt[3]
        l_sum_tab[ti]['baseline_5th'][stm] = bl_pcnt[1]
        l_sum_tab[ti]['baseline_80th'][stm] = bl_pcnt[2]

    # ---=== Get calibrators for the summary tables ===---
    # Initialise strings for each calibrator entry
    f_cal = ''
    p_cal = ''
    b_cal = ''
    # Iterate over each MS, add calibrator names to string ("|" delimited)
    for mspath in g_ms:
        # Add delimiter between each MS entry
        if f_cal != '':
            f_cal += '|'
        if p_cal != '':
            p_cal += '|'
        if b_cal != '':
            b_cal += '|'
        # Get calibrators
        cal_dict =  au.getCalibrators(vis=mspath,
                                      intent=['BANDPASS', 'FLUX', 'PHASE'],
                                      returnDict=True)
        # Add to strings
        f_cal += cal_dict['FLUX']
        p_cal += cal_dict['PHASE']
        b_cal += cal_dict['BANDPASS']
    # Add to summary tables
    c_sum_tab['flux_calibrator'][stm] = f_cal
    c_sum_tab['phase_calibrator'][stm] = p_cal
    c_sum_tab['bandpass_calibrator'][stm] = b_cal
    for ti in range(len(l_sum_tab)):
        l_sum_tab[ti]['flux_calibrator'][stm] = f_cal
        l_sum_tab[ti]['phase_calibrator'][stm] = p_cal
        l_sum_tab[ti]['bandpass_calibrator'][stm] = b_cal

    # ---=== Preparations for imaging continuum ===---
    band = 'continuum'
    '''
    NOTE: Continuum needs to be imaged first, as the continuum flux needs to be
          measured to check for possible continuum contamination in the line
          bands. If there is continuum above some threshold, then imcontsub()
          will be triggered after each line band has been imaged.
    '''
    # Create base name for image
    cont_im_name = g_name + '_' + band + '_' + t_arr
    cont_im_path = ppath + dir_name + cont_im_name # Absolute path to image

    # Add bands used for continuum to the summary table
    cb_str = ''
    for cb in cont_bands:
        if cb_str != '':
            cb_str += '|' # Add delimiter for subsequent entries
        cb_str += cb
    c_sum_tab['cont_bands'][stm] = cb_str

    # Get continuum spectral windows for imaging
    cont_spw_list = [] # List of all spw IDs for continuum
    cont_spw_12m = [] # List of 12m spw IDs for continuum
    for j, mspath in enumerate(g_ms):
        # Get dictionary of spws IDs and their band labels
        t_dict = au.transitions(mspath, source=iau_name[i],
                                intent='OBSERVE_TARGET#ON_SOURCE')
        spws = list(t_dict) # Get list of all spw IDs in MS
        spw_str = '' # Initialised string of spws for each MS
        # Get spw IDs for bands used for continuum imaging
        for sw in spws:
            for cb in cont_bands:
                if transitions_dict[cb] in t_dict[sw]:
                    if spw_str != '':
                        spw_str += ',' # Add commas before subsequent spws
                    spw_str +=str(sw)
                    '''
                    NOTE: This is a bit convoluted, but it works.
                    '''
        # Add string of spws for each MS to list
        cont_spw_list.append(spw_str)
        if g_ms_flag[j] == '12m':
            cont_spw_12m.append(spw_str)

    # ---=== Estimate maximum resolvable scale for continuum ===---
    cont_mrs_list = [] # List of MRS values
    for j, mspath in enumerate(g_ms):
        cont_mrs_list.append(MRS_L5(vis=mspath,
                                    field=iau_name[i],
                                    spw=cont_spw_list[j]))

    # Choose the largest MRS value of all MSs
    cont_mrs_max = np.max(cont_mrs_list) # This is in arcseconds

    # Add MRS to continuum summary table
    c_sum_tab['MRS_L5'][stm] = cont_mrs_max * u.arcsec

    # Set max tclean() component scales
    if limit_scales:
        l_scale = int(np.floor(cont_mrs_max / c_size))
    else:
        l_scale = int(-1)

    # ---=== Create dirty continuum image ===---
    # Initialise dictionary of clean parameters for dirty image
    dkwargs = dict()
    # Set tclean() parameters
    dkwargs.update(dict(vis=g_ms,
                        imagename=cont_im_path + '_dirty',
                        field=iau_name[i],
                        spw=cont_spw_list,
                        intent='OBSERVE_TARGET#ON_SOURCE',
                        specmode='mfs',
                        nterms=1,
                        gridder='mosaic',
                        mosweight=False,
                        imsize=im_size,
                        phasecenter=phase_center,
                        cell=cell_size,
                        deconvolver='asp',
                        largestscale=l_scale,
                        weighting=c_weighting,
                        robust=robust,
                        threshold='0.0mJy',
                        pbcor=False,
                        restoringbeam='common',
                        niter=0,
                        interactive=False
                        )
                   )

    # Run tclean() to create dirty continuum image
    tclean(**dkwargs)

    # Clean up unused by-products
    dirt_prods = glob.glob(ppath + dir_name + '*_' + band + '_*_dirty*')
    dirt_prods.remove(cont_im_path + '_dirty' + '.image')
    # Remove all other by-products of the dirty image creation
    for dirt_file in dirt_prods:
        shutil.rmtree(dirt_file)

    # ---=== Estimate noise threshold using the dirty image ===---
    # Clip pixels outside range [Q1 - fence*IQR, Q3 + fence*IQR]
    imstat_dict = imstat(cont_im_path + '_dirty' + '.image',
                         algorithm='hinges-fences', fence=1.5)

    calc_rms = imstat_dict['rms'][0] # Get RMS calculated by imstat()
    n_thresh = rms_factor * calc_rms # Multiply by the rms_factor
    threshold = str(n_thresh) + 'Jy' # Convert to string with units

    # ---=== Run tclean() to image continuum ===---
    # Prepare dictionary of tclean parameters
    kwargs = dict()
    # Set tclean() parameters
    kwargs.update(dict(vis=g_ms,
                       imagename=cont_im_path,
                       field=iau_name[i],
                       spw=cont_spw_list,
                       intent='OBSERVE_TARGET#ON_SOURCE',
                       specmode='mfs',
                       nterms=1,
                       gridder='mosaic',
                       mosweight=False,
                       imsize=im_size,
                       phasecenter=phase_center,
                       cell=cell_size,
                       deconvolver='asp',
                       largestscale=l_scale,
                       weighting=c_weighting,
                       robust=robust,
                       threshold=threshold,
                       pbcor=True,
                       restoringbeam='common',
                       niter=n_iter,
                       interactive=False,
                       fastnoise=False
                       )
                  )
    # Add automasking parameters
    kwargs.update(dict(automask_params[t_arr]['continuum']))

    # Run tclean with input parameter dictionary 'kwargs'
    tclean_dict = tclean(**kwargs)

    # Add flag for whether source was imaged
    c_sum_tab['imaged_flag'][stm] = True

    # Make tclean_dict more readily writable
    if isinstance(tclean_dict['summarymajor'], np.ndarray):
        tclean_dict['summarymajor'] = tclean_dict['summarymajor'].tolist()
    # Write dictionaries to file
    with open(cont_im_path + '_tclean_dict.json', 'w') as file:
       file.write(json.dumps(tclean_dict)) # Dictionary returned by tclean()
       file.close()
    with open(cont_im_path + '_tclean_kwargs.json', 'w') as file:
       file.write(json.dumps(kwargs)) # tclean() input arguments
       file.close()
    '''
    NOTE: These dictionaries are written out mostly for debugging purposes.
          See section "Returned Dictionary" for tclean() summary dictionary:
          https://casadocs.readthedocs.io/en/v6.5.4/
          notebooks/synthesis_imaging.html#Returned-Dictionary

          One thing that the tclean_kwargs file could be used for is to check if
          the current imaging parameters match exactly with imaging products
          that are already in the directory, and to skip re-cleaning if that is
          the case.
    '''
    # ---=== Add cleaned image properties to summary table ===---
    im_h = imhead(cont_im_path + '.image') # Header
    # Add beam properties
    b_p = im_h['restoringbeam']
    c_sum_tab['beam_major'][stm] = (b_p['major']['value'] *
                                     u.Unit(b_p['major']['unit']))
    c_sum_tab['beam_minor'][stm] = (b_p['minor']['value'] *
                                     u.Unit(b_p['minor']['unit']))
    c_sum_tab['beam_PA'][stm] = (b_p['positionangle']['value'] *
                                  u.Unit(b_p['positionangle']['unit']))
    # Get statistics
    im_s = imstat(cont_im_path + '.image') # Statistics
    im_sc = imstat(cont_im_path + '.image',
                   algorithm='hinges-fences', fence=1.5) # Statistics (clipped)
    # Add measured RMS
    im_rms = im_sc['rms'] * u.Unit(im_h['unit'])
    c_sum_tab['RMS'][stm] = im_rms
    # Add peak SNR and low SNR flag
    im_peak = im_s['max'] * u.Unit(im_h['unit'])
    c_sum_tab['peak_SNR'][stm] = im_peak / im_rms
    c_sum_tab['low_SNR_flag'][stm] = im_peak / im_rms < 5

    # Check for empty tclean() model
    md_s = imstat(cont_im_path + '.model') # Model statistics
    em_f = (md_s['max'] == 0) and (md_s['min'] == 0) and (md_s['sum'] == 0)
    c_sum_tab['empty_model_flag'][stm] = em_f

    # ---=== Export tlcean() continuum products as FITS ===---
    if export_fits:
        for ex in export_list:
            exportfits(imagename=cont_im_path + ex,
                       fitsimage=cont_im_path + ex + '.fits',
                       velocity=True, dropstokes=True, dropdeg=True,
                       optical=True, overwrite=True)

    # ---=== Image all specified line bands ===---
    # Iterate over each line band
    for li, l_band in enumerate(line_bands):
        # ---=== Preparations for imaging lines ===---
        band = l_band
        # Create base name for line image
        line_im_name = (g_name + '_' + band + '_' + str(ch_width) + 'kmps' +
                        '_' + t_arr)
        line_im_path = ppath + dir_name + line_im_name # Absolute path to image

        # Get current line spectral windows for imaging
        line_spw_list = [] # List of all spw IDs for current line band
        line_spw_12m = [] # List of 12m spw IDs for current line band
        for j, mspath in enumerate(g_ms):
            # Get dictionary of spws IDs and their band labels
            t_dict = au.transitions(mspath, source=iau_name[i],
                                    intent='OBSERVE_TARGET#ON_SOURCE')
            spws = list(t_dict) # Get list of all spw IDs in MS
            spw_str = '' # Initialised string of spws for each MS
            # Get spw IDs for bands used for selected line to be imaged
            for sw in spws:
                if transitions_dict[band] in t_dict[sw]:
                    if spw_str != '':
                        spw_str += ',' # Add commas before subsequent spws
                    spw_str +=str(sw)
            # Add string of spws for each MS to list
            line_spw_list.append(spw_str)
            if g_ms_flag[j] == '12m':
                line_spw_12m.append(spw_str)

        # ---=== Estimate maximum resolvable scale for current line band ===---
        line_mrs_list = [] # List of MRS values
        for j, mspath in enumerate(g_ms):
            line_mrs_list.append(MRS_L5(vis=mspath,
                                        field=iau_name[i],
                                        spw=line_spw_list[j]))

        # Choose the largest MRS value of all MSs
        line_mrs_max = np.max(line_mrs_list) # This is in arcseconds

        # Add MRS to continuum summary table
        l_sum_tab[li]['MRS_L5'][stm] = line_mrs_max * u.arcsec

        # Set max tclean() component scales
        if limit_scales:
            l_scale = int(np.floor(line_mrs_max / c_size))
        else:
            l_scale = int(-1)

        # ---=== Set channel width ===---
        width = str(ch_width) + 'km/s' # Convert float input to string

        # ---=== Create dirty image of current line band ===---
        # Initialise dictionary of clean parameters for dirty image
        dkwargs = dict()
        # Set tclean() parameters
        dkwargs.update(dict(vis=g_ms,
                            imagename=line_im_path + '_dirty',
                            field=iau_name[i],
                            spw=line_spw_list,
                            intent='OBSERVE_TARGET#ON_SOURCE',
                            specmode='cube',
                            width=width,
                            nterms=1,
                            gridder='mosaic',
                            mosweight=False,
                            imsize=im_size,
                            phasecenter=phase_center,
                            cell=cell_size,
                            deconvolver='asp',
                            largestscale=l_scale,
                            weighting=l_weighting,
                            robust=robust,
                            threshold='0.0mJy',
                            pbcor=False,
                            restoringbeam='common',
                            niter=0,
                            interactive=False
                            )
                       )

        # Run tclean() to create dirty continuum image
        tclean(**dkwargs)

        # Clean up unused by-products
        dirt_prods = glob.glob(ppath + dir_name + '*_' + band + '_' +
                               str(ch_width) + 'kmps_' + '*_dirty*')
        dirt_prods.remove(line_im_path + '_dirty' + '.image')
        # Remove all other by-products of the dirty image creation
        for dirt_file in dirt_prods:
            shutil.rmtree(dirt_file)

        # ---=== Estimate noise threshold using the dirty image ===---
        # Get total number of channels in dirty image
        dN_chan = imhead(line_im_path + '_dirty' + '.image')['shape'][-1]
        if l_Nch == 0:
            # Use last h_Nch channels of dirty cube
            chan_str = str(dN_chan - h_Nch) + '~' + str(dN_chan - 1)
        elif h_Nch == 0:
            # Use first l_Nch channels of dirty cube
            chan_str = '0~' + str(l_Nch)
        else:
            # Use first l_Nch and last h_Nch channels of dirty cube
            chan_str = ('0~' + str(l_Nch) + ',' +
                        str(dN_chan - h_Nch) + '~' + str(dN_chan - 1))
        imstat_dict = imstat(line_im_path + '_dirty' + '.image', chans=chan_str)

        calc_rms = imstat_dict['rms'][0] # Get RMS calculated by imstat()
        n_thresh = rms_factor * calc_rms # Multiply by the rms_factor
        threshold = str(n_thresh) + 'Jy' # Convert to string with units

        # ---=== Run tclean() to image current line band ===---
        # Prepare dictionary of tclean parameters
        kwargs = dict()
        # Set tclean() parameters
        kwargs.update(dict(vis=g_ms,
                           imagename=line_im_path,
                           field=iau_name[i],
                           spw=line_spw_list,
                           intent='OBSERVE_TARGET#ON_SOURCE',
                           specmode='cube',
                           width=width,
                           nterms=1,
                           gridder='mosaic',
                           mosweight=False,
                           imsize=im_size,
                           phasecenter=phase_center,
                           cell=cell_size,
                           deconvolver='asp',
                           largestscale=l_scale,
                           weighting=l_weighting,
                           robust=robust,
                           threshold=threshold,
                           pbcor=True,
                           restoringbeam='common',
                           niter=n_iter,
                           interactive=False,
                           fastnoise=False
                           )
                      )
        # Add automasking parameters
        kwargs.update(dict(automask_params[t_arr]['line']))

        # Run tclean with input parameter dictionary 'kwargs'
        tclean_dict = tclean(**kwargs)

        # Add flag for whether source was imaged
        l_sum_tab[li]['imaged_flag'][stm] = True

        # Make tclean_dict more readily writable
        if isinstance(tclean_dict['summarymajor'], np.ndarray):
            tclean_dict['summarymajor'] = tclean_dict['summarymajor'].tolist()
        # Write dictionaries to file
        with open(line_im_path + '_tclean_dict.json', 'w') as file:
           file.write(json.dumps(tclean_dict)) # Dictionary returned by tclean()
           file.close()
        with open(line_im_path + '_tclean_kwargs.json', 'w') as file:
           file.write(json.dumps(kwargs)) # tclean() input arguments
           file.close()
        '''
        NOTE: See above note about tclean() dictionaries below the continuum
              imaging tclean() call.
        '''

        # ---=== Continuum subtraction ===---
        # Measure peak continuum in flat image
        cont_stats = imstat(cont_im_path + '.image')
        cont_max = cont_stats['max']
        # Measure RMS in line-free channels of cube
        line_stats = imstat(line_im_path + '.image', chans=chan_str)
        line_rms = line_stats['rms']
        # Check continuum brightness against line RMS
        N_rms = 1 # N*RMS do you need continuum peak to be to trigger contsub()
        if cont_max >= N_rms * line_rms:
            run_contsub = True # Set flag to run imcontsub()
        else:
            run_contsub = False

        # Perform imconstub() if required
        if do_contsub & run_contsub:
            # Delete previous any previous continuum subtraction products
            cs_prods = glob.glob(line_im_path + '.cont*')
            for csp in cs_prods:
                if os.path.isfile(csp):
                    os.remove(csp)
                elif os.path.isdir(csp):
                    shutil.rmtree(csp)
            # Model and subtract continuum
            imcontsub(imagename=line_im_path + '.image',
                      linefile=line_im_path + '.contsub' + '.image',
                      contfile=line_im_path + '.contmod' + '.image',
                      fitorder=0,
                      chans=chan_str)
            '''
            NOTE: Consider limiting imcontsub() to a subregion where there is
                  only strong continuum. More channels could probably be used
                  to estimate the continuum also, so a method to reliably find
                  line-free channels automatically could be added earlier to the
                  pipeline. Ideally, this could be done on the dirty image as
                  then we could use more channels to estimate the noise
                  threshold for tclean(), but that's probably a bit
                  perfectionistic and will only be of marginal benefit.
            '''
            # Primary beam correct the continuum-subtracted image
            impbcor(imagename=line_im_path + '.contsub' + '.image',
                    pbimage=line_im_path + '.pb',
                    outfile=line_im_path + '.contsub' + '.image.pbcor',
                    overwrite=True)
            # Set flag in summary table
            l_sum_tab[li]['cont_sub_flag'][stm] = True

        # ---=== Add cleaned image properties to summary table ===---
        # Choose whether to get stats from normal or continuum subtracted image
        if do_contsub & run_contsub:
            file_ext = '.contsub.image'
        else:
            file_ext = '.image'
        im_h = imhead(line_im_path + file_ext) # Header
        # Add beam properties
        b_p = im_h['restoringbeam']
        l_sum_tab[li]['beam_major'][stm] = (b_p['major']['value'] *
                                            u.Unit(b_p['major']['unit']))
        l_sum_tab[li]['beam_minor'][stm] = (b_p['minor']['value'] *
                                            u.Unit(b_p['minor']['unit']))
        l_sum_tab[li]['beam_PA'][stm] = (b_p['positionangle']['value'] *
                                         u.Unit(b_p['positionangle']['unit']))
        # Get statistics
        im_s = imstat(line_im_path + file_ext) # Statistics
        im_slf = imstat(line_im_path + file_ext, chans=chan_str) # Line-free
        # Add measured RMS
        im_rms = im_slf['rms'] * u.Unit(im_h['unit'])
        l_sum_tab[li]['RMS'][stm] = im_rms
        # Add peak SNR and low SNR flag
        im_peak = im_s['max'] * u.Unit(im_h['unit'])
        l_sum_tab[li]['peak_SNR'][stm] = im_peak / im_rms
        l_sum_tab[li]['low_SNR_flag'][stm] = im_peak / im_rms < 5

        # Check for empty tclean() model
        md_s = imstat(line_im_path + '.model') # Model statistics
        em_f = (md_s['max'] == 0) and (md_s['min'] == 0) and (md_s['sum'] == 0)
        l_sum_tab[li]['empty_model_flag'][stm] = em_f

        # ---=== Export tlcean() cube products as FITS ===---
        if export_fits:
            if do_contsub & run_contsub:
                # Use different export list for continuum-subtracted cubes
                cs_ex_list = (set(export_list) -
                              set(['.image', '.image.pbcor']))
                cs_ex_list = list(cs_ex_list |
                              set(['.contsub.image', '.contsub.image.pbcor']))
                for ex in cs_ex_list:
                    exportfits(imagename=line_im_path + ex,
                               fitsimage=line_im_path + ex + '.fits',
                               velocity=True, dropstokes=True, dropdeg=True,
                               optical=True, overwrite=True)
            else:
                for ex in export_list:
                    exportfits(imagename=line_im_path + ex,
                               fitsimage=line_im_path + ex + '.fits',
                               velocity=True, dropstokes=True, dropdeg=True,
                               optical=True, overwrite=True)
        '''
        NOTE: I don't like how contrived this is, but it's what I could do in
              the time that I had.
        '''

        # ---=== Convert FITS cubes from Jy/beam -> K ===---
        # Only do this if export_fits = True
        if export_fits:
            # List of file types to convert
            if do_contsub & run_contsub:
                conv_list = ['.contsub.image', '.contsub.image.pbcor']
            else:
                conv_list = ['.image', '.image.pbcor']
            # Import each FITS file using astropy
            for conv in conv_list:
                fits_d = fits.getdata(filename=line_im_path + conv + '.fits',
                                      ext=0)
                fits_h = fits.getheader(filename=line_im_path + conv + '.fits',
                                      ext=0)
                # Check unit and do conversion
                if u.Unit(fits_h['BUNIT']) == (u.Jy/u.beam):
                    # Calculate beam area
                    omega_B = ((np.pi / (4 * np.log(2))) *
                               fits_h['BMAJ'] * u.deg * fits_h['BMIN'] * u.deg)
                    # Get rest frequency of data
                    rest_nu = fits_h['RESTFRQ'] * u.Hz
                    fits_d = fits_d << u.Unit(fits_h['BUNIT']) # Apply unit
                    # Define equivalency
                    eqv = u.brightness_temperature(frequency=rest_nu,
                                                   beam_area=omega_B)
                    # Convert to brightness temperature K
                    fits_d = fits_d.to(u.K, equivalencies=eqv)
                    # Update header, strip unit from data
                    fits_h['BUNIT'] = fits_d.unit.to_string()
                    fits_d = fits_d.value
                    # Write FITS file with updated units
                    fits.writeto(filename=line_im_path + conv + '.fits',
                                 data=fits_d, header=fits_h, overwrite=True)
        '''
        NOTE: This is not very read/write efficient, but it is a lot easier
              than doing the conversion in CASA before exporting the cubes as
              FITS files. This could be merged with the above block of code and
              cleaned up as it is equally contrived.
        '''

# ---=== Write out summary tables ===---
# Write out continuum summary table
c_sum_tab.write(ppath + dname + '/original/' + cont_tab_name,
                format='fits', overwrite=True)
# Read in table again to add TOPCAT-readable column descriptions
ftab = fits.open(ppath + dname + '/original/' + cont_tab_name)
# Get header
header = ftab[1].header
# Iterate over each header key, add column description to TCOMM key
for key in header.keys():
    if key.startswith('TTYPE'):
        key_id = key.strip('TTYPE')
        key_str = header[key]
        header['TCOMM' + key_id] = c_sum_tab[key_str].description
# Write to file
ftab.writeto(ppath + dname + '/original/' + cont_tab_name,
             overwrite=True)

# Write out line summary tables
for ti, lst in enumerate(l_sum_tab):
    lst.write(ppath + dname + '/original/' + line_tab_names[ti],
              format='fits', overwrite=True)
    # Read in table again to add TOPCAT-readable column descriptions
    ftab = fits.open(ppath + dname + '/original/' + line_tab_names[ti])
    # Get header
    header = ftab[1].header
    # Iterate over each header key, add column description to TCOMM key
    for key in header.keys():
        if key.startswith('TTYPE'):
            key_id = key.strip('TTYPE')
            key_str = header[key]
            header['TCOMM' + key_id] = lst[key_str].description
    # Write to file
    ftab.writeto(ppath + dname + '/original/' + line_tab_names[ti],
                 overwrite=True)
'''
NOTE: Some FITS table viewing applications (e.g. TOPCAT) aren't able to show the
      column descriptions that astropy adds to the FITS files, so this is a
      workaround so that the descriptions show up in TOPCAT.
'''
