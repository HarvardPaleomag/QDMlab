Release Notes
*************
2021.1.4
========
  features
  --------
  - adds option alignmentType: str ['laser']
    allows `subtract_blank` to subtract any two maps
    use `alignmentType, 'led'`
  - adds polynomial option for image alignment
  - fixes issue with misaligned maps after fitting

  important changes
  -----------------
  - function `dipole_fit` renamed now `fit_source`
  - function `dipole_fit_series` renames now `fit_sources_series`

  minor fixes
  -----------
  - removes use of `nanstd`, `nanmean` and `nanmedian` replaced with `std`, `mean`,`median`
  - small fixes and optimization in `ODMR_to_B111_plot`, `QDM_figure`, `scalebar` and `demag_behavior_plot`
  - removes `reshape_QDM_data`

2021.1.3
========
  features
  --------
  - adds python admin_tools
  - adds `projective`, `polynomial` options to `get_image_tform_complex`
  - adds colormap switching to `QDM_figure`
  - adds `ODMR_to_B111_plot`
  - adds more information to return value of `dipole_fit_series`
  - adds title option to `crop_map`
  - adds option not to threshold in `filter_hot_pixels`

  important changes
  -----------------
  - `demag_behavior` keyword `pixelError` now called `pixelShift`
  - ODMR_to_B111 now returns all fits

  minor fixes
  -----------
  - fixes issue with `tform_data` if binning wasn't correct
  - adds frequencies xlabels in globalFraction_estimator
  - proper logging in `subtract_constant`
  - indexing in `xy2index` and `index2xy`

2021.1.2a hotfix
================
- fixes bug where an error would occur in `QDM_figure` if all values are 0 or nan
- fixes typo in `get_transformed_maps` when using `'checkPlot'` keyword
- fixes typo in `pick_box`
- fixes bug in `show_references`

2021.1.2
========
  features
  --------
  - `dipole_fit_series` now can use the constrained code by setting
    `'constrained', true` and passing `m0`, `hguess`, `minheight`,
    `maxheight`, `boxwidth` to the desired values

  important changes
  -----------------
  - `QDM_lorentzian_fit` is *deprecated* now. Use `ODMR_to_B111`
  - `estimate_coercivity` is *deprecated* now. Use `demag_behavior`
  - `coercivity_results_plot` is *deprecated* now. Use `demag_behavior_plot`

  minor fixes
  -----------
  - fixes issue where QDM_figure does not show the data
  - fixes issue where maps are displayed in only two colors
  - slopeCorrection in `QDMR_to_B111` now is global to all types and affects the GPUdata

2021.1.1
========
  - adds precompiled GPUfit package with corresponding functions.
    **NOTE:** ESR3RT was renamed to ESRN14
  - adds GPUfit modified files for compiling
  - adds `scalebar` function to add a scale-bar to the map (QDM_figure)
  - adds warning pop-up if NaN values are detected in B111 data, when converting to Bz
  - filterProps can now be passed as a structure into `estimate_coercivity` and `QDM_figure`

  minor fixes
  -----------
  - fixes bug where binning can not be automatically determined
  - adds title to viscosity plots
  - better naming of resulting maps from `subtract_source`
  - fixes issue where LED would be cropped after subtract_blank
  - fixes issue with wrong cropping of LED in QDMdataprocessing
  - `coercivity_results_plot` now has option to add mean for all ROIs
  - better first line comments
  - replaces `exist_struct` with `isfield`

2021.0.beta4
============
- changes default globalFraction to 0.25
- adds tracking of failed pixels
- QDM_lorentian_fits now saves all the fitting related data (i.e. the return value) to 'final_fits_(binxbin).mat'
