<p align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="docs/Logo_dark.png">
        <img alt="Shows a logo with black text for light colour themes and white text in dark colour themes." src="docs/Logo_light.png">
    </picture>
</p>

# KILOGAS Imaging
This repo contains all scripts relating to processing ALMA data for the KILOGAS project.

## Requirements
These scripts were written and tested on Measurements Sets (MSs) from ALMA Cycle 11, using [CASA 6.5.4.9](https://almascience.nrao.edu/processing/science-pipeline), and make extensive use of functions from [Analysis Utilities](https://safe.nrao.edu/wiki/bin/view/Main/CasaExtensions) by [Hunter et al. (2023)](https://zenodo.org/records/13887809).

Please follow the recommended installation instructions for [CASA](https://casadocs.readthedocs.io/en/stable/notebooks/introduction.html#id1) and [Analysis Utilities](https://casaguides.nrao.edu/index.php/Analysis_Utilities).

The scripts will likely work with other later (and possibly some earlier) versions of CASA and MSs from different ALMA Cycles, but this has not been tested.

## Setup
The scripts require calibrated MSs that are located within the standard directory structure that is generated when running the pipeline calibration from the ALMA Science Archive, i.e. each calibrated MS will be located in:
```
[project code]/science_goal.uid___*/group.uid___*/member.uid___*/calibrated/
```
This was a design choice motivated by not making the user have to manage the location of MSs and perform lots of spliting and concatenating prior to imaging, as these both take a loot of time and can increase the size on disk substantially (especially for large projects). It also makes it easier to trace any issues in the imaging products back to the source MSs.

### band_names.json
This is a dictionary containing the long band names and their respective short names (user defined) that are used to refer to the bands in the rest of the scripts and iamging products. The long names can be found by running
```python
CASA <1>: au.transitions(vis='/abs/path/to/example.ms')
```
which will be consistent across the entire ALMA project.

### config.ini
This is the main confguration file that all of the scripts depend on. Default values are shown in square brackets. It is broken up into four main sections:
#### paths
* `data_path` - Path to top-level directory of ALMA data (usually the project code).
* `prod_path` - Path to top-level directory of imaging products.
* `dir_name` - Name to give top-level directory of output (e.g. DR1, test_01).

#### input
* `targets` - Target selection, can be either `all` (which will select all targets that have been found within the project), a comma-delimited list of integers corresponding to their KILOGAS IDs, or a the name of `.txt` file containing a comma-delimited list of KILOGAS IDs. `[all]`
* `line_bands` - Which bands are to be imaged as line emission datacubes, which can be a list of multiple bands. These need to be consistent with the short names set in `band_names.json`. `[co2-1]`
* `cont_bands` - Which bands are to be jointly imaged using multifrequency synthesis (mfs) for the continuum. These need to be consistent with the short names set in `band_names.json`. `[continuum,cs5-4,ch3oh3-2]`

#### imaging
* `cell_size` - Image cell size in arcseconds. `[0.1]`
* `rms_threshold` - Noise threshold for stopping `tclean()`, in multiples of the measured RMS. `[1]`
* `chan_width` - Channel width for line imaging in km/s. `[10.0]`
* `nchan_rms` - Number of channels from the start/end of the spectrum to use for RMS measurement when line imaging. Can be a single value, or two values corresponding the first N channels and last M channels of the cube. `[10]`
* `line_weight` - Weighting scheme for line imaging. `[briggs]`
* `cont_weight` - Weighting scheme for continuum imaging. `[natural]`
* `robust` - Robust parameter value for when Briggs weighting is used. `[0.5]`
* `limit_scales` - Limit the maximum CLEAN component scale used by the ASP algorithm. This will be set to the estimated maximum resolvable scale (MRS_L5). `[False]`
* `do_contsub` - Allow the imaging pipeline to perform continuum subtraction of the line cubes where necessary. `[True]`

#### output
* `export_fits` - Export FITS files of imaging products? `[True]`
* `fits_list` - List of image types to export as FITS. If `default` then use predefined list: .image, .image.pbcor, .residual, .pb, .model, .mask. `[default]`
* `cube_in_K` - Convert line cubes from Jy/beam to K? `[True]`
* `use_optical` - Use optical velocity convention in FITS files? `[True]`

# Usage
There are three main scripts that are all are intended to be run within CASA.

## ms_finding.py
This script is responsible for locating all of the paths to the MSs in the project and matching them to the science targets withing them. This script only needs to be run once, or if more targets/observations are added to the project or any changes are made to the paths of the MSs. It requires a FITS table that contains the KILOGAS IDs, IAU names and ALMA Scheduling Block numbers for all targets (this is located in `tables/KILOGAS_global_catalog.fits`). To run this script, start CASA in the parent directory of this repo and then type:
```python
CASA <1>: execfile('kilogas_imaging/ms_finding.py')
```

It will list any targets it cannot find MSs for, and produce a FITS table (locate in `tables/kilogas_ms_paths.fits`) that contains:
* `KGAS_ID` - KILOGAS ID
* `IAUname` - IAU name
* `SB_num` - Scheduling Block (SB) number
* `has12m?` - Boolean flag for if 12m MS is present
* `12m_path` - List of 12m MS paths, "|" delimited
* `has7m?` - Boolean flag for if 7m MS is present
* `7m_path` - List of 7m MS paths, "|" delimited


## kilogas_imaging.py
This is the pipeline responsible for imaging the ALMA data. It automatically calculated and sets the vast majority of parameters required for `tclean()`, and only depends of the parameters set in `config.ini`. The pipeline needs `tables/kilogas_ms_paths.fits` to run.

To run the imaging pipeline, start CASA in the parent directory of this repo and then type:
```python
CASA <2>: execfile('kilogas_imaging/kilogas_imaging.py')
```

Here is a brief summary of the pipeline's operation:
1. Read in and parse the contents of `config.ini`, `band_names.json`, `automask_params.json`, and `tables/kilogas_ms_paths.fits`
2. Iterate over each of the sources selected for imaging:
  - Check for any conflicting imaging products in the output directory, and remove them.
  - Calculate phasecentre, optimal image size
  - Continuum imaging:
    - Get spectral windows for continuum imaging
    - Generate a dirty image of the continuum
    - Estimate the noise threshold for `tclean()` using the dirty image, by calculating the RMS using the Tukey's fences method (i.e. `algorithm=hinges-fences` in `imstat()`)
    - Perform continuum imaging
    - Export continuum imaging products, and JSON files of `tclean()` input/output (useful for debugging)
  - Line imaging (i.e. iterate over each line band in the config):
    - Get spectral windows for line imaging
    - Generate a dirty image cube of the line
    - Estimate the noise threshold for `tclean()` using the dirty image, using the RMS in the channels specified in the config
    - Perform line imaging
    - Check if measured continuum is >1*RMS in the imaged line cube
      - If yes, perform continuum subtraction using `imcontsub()`
    - Convert line image cubes from units of Jy/beam to K
    - Export line imaging products, and JSON files of `tclean()` input/output (useful for debugging)
3. Write out summary tables for each imaged band, with a row for each imaged source.

The imaging product names follow the format:
```
KGAS[ID]_[band]{_channel width}_[telescope array].[image type]
```
They are written to `[prod_path]/[dir_name]/original/KGAS[ID]/`, and the summary tables are written to `[prod_path]/[dir_name]/`. Note that the channel width is only included for line image cubes, and not for the continuum images.

The summary tables produced by the pipeline contain quite a lot of useful information, such as:
* Beam properties (major and minor axis lengths, position angle)
* Measured RMS
* Peak SNR
* Maximum resolvable scale
* Projected baseline lengths (min, max, 5th and 80th percentile)
* Calibrators (flux, phase, bandpass)
If the pipeline is re-run then these tables will be read in and updated, so if only you want to re-image only a few sources, only those sources will be updated in the summary tables.

## ifu_matching.py
This script is responsible for convolving and regridding the ALMA line products to the resolution and grid size of the MaNGA/SAMI IFU data. This script requires a table (`ifu_info/matched_ids_fwhm.txt`) that contains KILOGAS IDs, MaNGA/SAMI ID's and the FWHM of the IFU data for each source. The script also required some a single-plane reference IFU image for the regridding step, whicha are located in `ifu_info/ifu_data/`.

To run this script, start CASA in the parent directory of this repo and then type:
```python
CASA <3>: execfile('kilogas_imaging/ifu_matching.py')
```

This script will iterate over all sources specified in `config.ini`, and then:
1. Import the FITS files for the non-PB corrected ALMA line cube, the ALMA primary beam, and the corresponding IFU reference image into CASA image format, and get the IFU FWHM for that field
2. Smooth the ALMA line cube to the IFU target FWHM
3. Primary beam correct the smoothed ALMA cube
4. Regrid the smoothed and PB-corrected smoothed ALMA line cubes to the IFU grid
5. Export the IFU-matched files as FITS (appended with `.ifumatched.fits`)

The IFU-matched products are written to `[prod_path]/[dir_name]/matched/KGAS[ID]/`.
