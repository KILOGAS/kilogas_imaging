"""
Python script that locates sources within their respective measurement sets,
maps their location within the project's directory structure, and outputs a
FITS table with rows for each KILOGAS source.
Tested on CASA 6.5.4 onwards and MS from ALMA Cycle 11.
"""

import configparser
import glob
import os
import sys

import analysisUtils as au
import numpy as np

from astropy.table import Table, Column
from astropy.io import fits
# from tqdm import tqdm

# ==============================================================================

# ---=== Define paths, import catalog ===---
# Get location of repository
rpath = os.path.abspath(__file__).strip('ms_finding.py')

# Get other paths from the main config file
cfg = configparser.ConfigParser(inline_comment_prefixes='#')
cfg.read(rpath + 'config.ini')
# Absolute path to ALMA data folder
dpath = cfg.get('paths', 'data_path')

# Import global catalog of KILOGAS sources
kgas_cat = Table.read(rpath + 'tables/KILOGAS_global_catalog.fits',
                      format='fits')

# ---=== Get paths of all MSs ===---
# Locate all MSs in top-level directory, and create list of all paths
glob_base = '*/*/*/' + 'calibrated/' # Base directory pattern to search
glob_ms = glob_base + '*.ms' # All MSs
glob_targets = glob_base + '*_targets.ms' # MSs of targets only
glob_line = glob_base + '*_targets_line.ms' # MSs of cont. subtracted targets

# Remove all targets.ms and targets_line.ms from list of all MSs
base_ms_set = (set(glob.glob(dpath + glob_ms)) -
               set(glob.glob(dpath + glob_targets)) -
               set(glob.glob(dpath + glob_line)))

# Convert to list of only un-split MSs, i.e. uid___####_########_######.ms
base_ms = list(base_ms_set)

# ---=== Initialise output table ===---
# Create table containing only relevant columns from
kgas_tab = kgas_cat['KGAS_ID', 'IAUname', 'SB_num']
# Convert MaskedColumn to Column
kgas_tab['SB_num'] = kgas_tab['SB_num'].filled()
# Convert to unicode strings make string comaprison easier later on
kgas_tab.convert_bytestring_to_unicode()

N_source = len(kgas_tab) # Get number of targets

# Initialise columns to be added to table (lists will be converted later)
has_12m = Column(np.zeros(N_source).astype(bool), name='has12m?') # 12m flag
path_12m = [''] * N_source # List of 12m MS paths
has_7m = Column(np.zeros(N_source).astype(bool), name='has7m?') # 7m flag
path_7m = [''] * N_source # List of 7m MS paths

# ---=== Populate table ===---
kg_set = set(kgas_tab['KGAS_ID']) # Set of all KGAS sources
ms_set = set() # Set to keep track of all KGAS sources with 12m MSs found

# Iterate over list of MSs
# ms_pbar = tqdm(total=len(kg_set), desc='Targets found in 12m MSs',
               # leave=False) # Progress bars for found MSs
# for i, m in enumerate(tqdm(base_ms, desc='MSs checked', leave=False)):
'''
NOTE: tqdm is not part of the base CASA installation, so this is commented out
      so that no additional installations are required by the user.
'''
for i, m in enumerate(base_ms):
    # Find science targets within the current MS
    targets = au.getTargetsForIntent(vis=m, intent='OBSERVE_TARGET#ON_SOURCE')
    # Get indices of found science targets in table
    targ_inds = np.where(np.isin(kgas_tab['IAUname'], np.array(targets)))[0]
    kg_ids = kgas_tab['KGAS_ID'][targ_inds] # Get relevant KGAS_IDs from table

    # Track KGAS_IDs found in the MSs
    old_len = len(ms_set) # Measure length of set prior to adding new KGAS_IDs

    # Get array config from MS
    arr_conf = au.getConfig(m)
    # Set flags for array type
    is7m = arr_conf == '7M'
    is12m = arr_conf != '7M'

    # Fill in flags and paths for 12m and 7m MSs
    if is7m:
        has_7m[targ_inds] = True
        for ti in targ_inds: # Iterate because of variable length strings
            # If path field already filled, if so append path with "|" delimiter
            if path_7m[ti] == '':
                path_7m[ti] = m
            else:
                path_7m[ti] += '|' + m
    if is12m:
        has_12m[targ_inds] = True
        for ti in targ_inds: # Iterate because of variable length strings
            # If path field already filled, if so append path with "|" delimiter
            if path_12m[ti] == '':
                path_12m[ti] = m
            else:
                path_12m[ti] += '|' + m
        # Add KGAS_IDs to set (if it is found in a 12m MS)
        ms_set.update(set(kg_ids))

    # Update MS progress bar (tqdm)
    # ms_pbar.update(len(ms_set) - old_len)
    # Update progress using native python
    sys.stdout.write('\rMSs checked: {}/{}. '.format(i+1, len(base_ms)) +
                     'Targets found in 12m MSs: ' +
                     '{}/{}'.format(len(ms_set), N_source))
    sys.stdout.flush()

# Print summary of 12m MSs found
print('\n\n12m MSs found for ' +
      '{0}/{1} KGAS targets.'.format(len(ms_set), len(kg_set)))
# Print list of all KGAS_IDs where MSs are not found
noms_set = kg_set - ms_set
if len(noms_set) > 0:
        print('\nNo 12m MSs found for ' +
              '{0} KGAS targets:\n{1}'.format(len(noms_set), sorted(noms_set)))

# ---=== Format table for output and write to file ===---

# Convert all lists to columns
p12m = Column(data=path_12m, name='12m_path')
p7m = Column(data=path_7m, name='7m_path')

# Append all new columns to kgas_tab
kgas_tab.add_columns(cols=[has_12m, p12m, has_7m, p7m])
kgas_tab.convert_unicode_to_bytestring() # Convert back to bytestring

# Add descriptions to table columns
kgas_tab['KGAS_ID'].description = 'KILOGAS ID'
kgas_tab['IAUname'].description = 'IAU name'
kgas_tab['SB_num'].description = 'Scheduling Block (SB) number'
kgas_tab['has12m?'].description = 'Flag for if 12m MS is present'
kgas_tab['12m_path'].description = 'List of 12m MS paths, "|" delimited'
kgas_tab['has7m?'].description = 'Flag for if 7m MS is present'
kgas_tab['7m_path'].description = 'List of 7m MS paths, "|" delimited'

# Write table to file
kgas_tab.write(rpath + 'tables/kilogas_ms_paths.fits',
               format='fits', overwrite=True)

# Read in table again to add TOPCAT-readable column descriptions
ftab = fits.open(rpath + 'tables/kilogas_ms_paths.fits')
# Get header
header = ftab[1].header
# Iterate over each header key, add column description to TCOMM key
for key in header.keys():
    if key.startswith('TTYPE'):
        key_id = key.strip('TTYPE')
        key_str = header[key]
        header['TCOMM' + key_id] = kgas_tab[key_str].description
# Write to file
ftab.writeto(rpath + 'tables/kilogas_ms_paths.fits', overwrite=True)
'''
NOTE: Some FITS table viewing applications (e.g. TOPCAT) aren't able to show the
      column descriptions that astropy adds to the FITS files, so this is a
      workaround so that the descriptions show up in TOPCAT.
'''
