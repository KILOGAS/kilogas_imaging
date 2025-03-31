"""
Script that regrids KILOGAS ALMA data products to match MaNGA/SAMI IFU data.
Tested on CASA 6.5.4.
"""

import configparser
import glob
import os
import re
import shutil
import sys

import astropy.units as u

from astropy.table import Table, Column
from astropy.io import fits

# ==============================================================================

# ---=== Define paths, import catalog ===---
# Get location of repository
rpath = os.path.abspath(__file__).strip('ifu_matching.py')

# Get other paths from the main config file
cfg = configparser.ConfigParser(inline_comment_prefixes='#')
cfg.read(rpath + 'config.ini')
# Absolute path to ALMA product folder
ppath = cfg.get('paths', 'prod_path')
# Name of imaging product output directory
dname = cfg.get('paths', 'dir_name').strip('/') + '/'
# Target selection
source_list = cfg.get('input', 'targets')
if source_list=='all':
    # Image all available targets with 12m MSs in table
    source_list = 'all'
else:
    # Convert string to list of targets ids (integers)
    source_list = [int(id) for id in re.split(', |,| ',
                                              source_list.strip('()[]{}'))]
# Bands to regrid
line_bands = re.split(', |,| ', cfg.get('input', 'line_bands').strip('()[]{}'))
# Channel width
ch_width = cfg.getfloat('imaging', 'chan_width')
'''
NOTE: This probably should be wrapped up in a function like the imaging script
'''

# Import table of KILOGAS IDs, IFU IDs and IFU FWHM values
ifu_tab = Table.read(rpath + 'ifu_info/matched_ids_fwhm.txt',
                     format='ascii')

# ---=== Prepare summary table ===---
# Survey name
srv_name = ifu_tab['Survey']
srv_name[srv_name=='MANGA'] = 'MaNGA' # Change to match filename
# Target ID numbers for their respective survey
srv_id = ifu_tab['ID']
srv_id.name = 'Survey_ID' # Rename for consistency
# KILOGAS ID
kgas_id = ifu_tab['KGAS_ID']
# IFU FWHM
ifu_fwhm = ifu_tab['FWHM']
ifu_fwhm = ifu_fwhm * u.arcsec # Add unit to column
# Create column that shows whether the data have been matched or not
r_flag = Column(name='regrid_flag', length=len(kgas_id), dtype=bool,
                description='True if IFU matching performed')
# Initialise summary table
ifu_sum_tab = Table()
ifu_sum_tab.add_columns([kgas_id, srv_name, srv_id, ifu_fwhm,r_flag])
# Add column descriptions
ifu_sum_tab['KGAS_ID'].description = 'KILOGAS ID'
ifu_sum_tab['Survey'].description = 'Survey name (MaNGA/SAMI)'
ifu_sum_tab['Survey_ID'].description = 'Target ID (MaNGA/SAMI)'
ifu_sum_tab['FWHM'].description = 'FWHM of the IFU data (MaNGA/SAMI)'
# Sort table by KILOGAS ID
ifu_sum_tab.sort('KGAS_ID')

# Create source list if it was set to "all" in the config
if source_list == 'all':
    source_list = ifu_sum_tab['KGAS_ID'].tolist()

# ---=== Run IFU matching procedure ===---
# Iterate over KGAS_IDs
for i, kg_id in enumerate(source_list):
    # Create string of KILOGAS ID
    g_name = 'KGAS' + str(kg_id)
    # Create directory structure
    dir_name = dname + 'matched/' + g_name + '/'
    # Check if file already exists
    if os.path.exists(ppath + dir_name):
        for lb in line_bands:
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
    else:
        # Create directory
        os.makedirs(ppath + dir_name)
    # Create mask for IFU table
    ifu_mask = ifu_sum_tab['KGAS_ID'] == kg_id
    # Import IFU data
    ifu_name = ('KG-' + ifu_sum_tab['Survey'][ifu_mask][0] + '-' +
                ifu_sum_tab['Survey_ID'][ifu_mask][0] +
                '.Pipe3D.cube.fits')
    ifu_file = rpath + 'ifu_info/ifu_data/' + ifu_name
    ifu_casa = ppath + dir_name + g_name + '_ifu.image'
    importfits(fitsimage=ifu_file, imagename=ifu_casa,
               zeroblanks=False, overwrite=True)
    # Get target beamsize for smoothing
    target_fwhm = str(ifu_sum_tab['FWHM'][ifu_mask][0]) + 'arcsec'
    # Iterate over each line band
    for lb in line_bands:
        # Find files to import
        base_path = ppath + dname + 'original/' + g_name + '/'
        # FITS version of *.image cube (K converted, contsub when needed)
        im_path = (base_path + g_name + '_' +  lb + '_' +
                   str(ch_width) + 'kmps_' + '*.image.fits')
        # FITS version of *.pb cube
        pb_path = (base_path + g_name + '_' +  lb + '_' +
                   str(ch_width) + 'kmps_' + '*.pb.fits')
        im_file = glob.glob(im_path)
        pb_file = glob.glob(pb_path)
        if len(im_file) > 0:
            # Import ALMA image file
            temp_im_file = (ppath + dir_name +
                            im_file[0].split('/')[-1].strip('.fits'))
            importfits(fitsimage=im_file[0], imagename=temp_im_file,
                       zeroblanks=False, overwrite=True)
            # Smooth ALMA image to target resolution
            imsmooth(imagename=temp_im_file, outfile=temp_im_file + '.smoothed',
                     kernel='gauss', major=target_fwhm, minor=target_fwhm,
                     pa='0.0deg', targetres=True, overwrite=True)
            # Import PB
            temp_pb_file = (ppath + dir_name +
                            pb_file[0].split('/')[-1].strip('.fits'))
            importfits(fitsimage=pb_file[0], imagename=temp_pb_file,
                       zeroblanks=False, overwrite=True)
            # Primary beam correct the smoothed image
            impbcor(imagename=temp_im_file + '.smoothed',
                    pbimage=temp_pb_file,
                    outfile=temp_im_file + '.smoothed' + '.pbcor',
                    overwrite=True)
            # Regrid smoothed image to IFU footprint
            ifu_match_name = temp_im_file + '.ifumatched'
            imregrid(imagename=temp_im_file + '.smoothed',
                     output=ifu_match_name,
                     axes=[0, 1], template=ifu_casa)
            # Export IFU matched image as FITS
            exportfits(imagename=ifu_match_name,
                       fitsimage=ifu_match_name + '.fits',
                       velocity=True, dropstokes=True, dropdeg=True,
                       optical=True, overwrite=True)
            # Regrid PB-corrected smoothed image to IFU footprint
            ifu_match_name = temp_im_file + '.pbcor' + '.ifumatched'
            imregrid(imagename=temp_im_file + '.smoothed' + '.pbcor',
                     output=ifu_match_name,
                     axes=[0, 1], template=ifu_casa)
            # Export IFU matched image as FITS
            exportfits(imagename=ifu_match_name,
                       fitsimage=ifu_match_name + '.fits',
                       velocity=True, dropstokes=True, dropdeg=True,
                       optical=True, overwrite=True)
            # Set regrid_flag to True in summary table
            ifu_sum_tab['regrid_flag'][ifu_mask] = True
        else:
            print('No image file found for ' + g_name)
    # Update progress using native python
    sys.stdout.write('\rTargets regridded: ' +
                     '{}/{}.'.format(i + 1, len(source_list)))
    sys.stdout.flush()

# Write summary table to file
ifu_sum_tab.write(ppath + dname + 'matched/' + 'ifu_matching_table.fits',
                  format='fits', overwrite=True)
'''
NOTE: The name of this table should be constructed a bit more carefully (like
      in the imaging script) so that you can have one per band/channel width
      that you run. That being said, this table is mostly for debugging the run you're currently doing, so not a super urgent fix.
'''

# Read in table again to add TOPCAT-readable column descriptions
ftab = fits.open(ppath + dname + 'matched/' + 'ifu_matching_table.fits')
# Get header
header = ftab[1].header
# Iterate over each header key, add column description to TCOMM key
for key in header.keys():
    if key.startswith('TTYPE'):
        key_id = key.strip('TTYPE')
        key_str = header[key]
        header['TCOMM' + key_id] = ifu_sum_tab[key_str].description
# Write to file
ftab.writeto(ppath + dname + 'matched/' + 'ifu_matching_table.fits',
             overwrite=True)
'''
NOTE: Some FITS table viewing applications (e.g. TOPCAT) aren't able to show the
      column descriptions that astropy adds to the FITS files, so this is a
      workaround so that the descriptions show up in TOPCAT.
'''
