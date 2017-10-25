# MATLAB/Fieldtrip Yokogawa Coregistration Scripts

## coreg_yokogawa_icp is a work in progress scripts to improve the coregistation between structural MRI and MEG data with polhemus headshape data. This approach uses the iterative closest point (ICP) algorithm to match scalp surface with downsampled headshape information.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coreg_yokogawa_icp is a function to coregister a structural MRI with MEG data 
% and associated headshape information
%
% Written by Robert Seymour Oct 2017 (some subfunctions written by Paul
% Sowman)
%
% INPUTS:
% - dir_name        = directory name for the output of your coreg
% - confile         = full path to the con file
% - mrkfile         = full path to the mrk file
% - mri_file        = full path to the NIFTI structural MRI file
% - hspfile         = full path to the hsp (polhemus headshape) file
% - elpfile         = full path to the elp file
% - hsp_points      = number of points for downsampling the headshape (try 100-200)
% - scalpthreshold  = threshold for scalp extraction (try 0.05 if unsure)
%
% OUTPUTS:
% - grad_trans              = correctly aligned sensor layout 
% - headshape_downsampled   = downsampled headshape (original variable name I know) 
% - mri_realigned           = the mri realigned based on fiducial points
% - trans_matrix            = transformation matrix for accurate coregistration
% - headmodel_singleshell   = coregistered singleshell headmodel
% 
% THIS IS A WORK IN PROGRESS FUNCTION - any updates or suggestions would be
% much appreciated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

![Imgur](https://i.imgur.com/UDnlAqA.png)

## Papers:

- [Gohel, B., Lim, S., Kim, M. Y., Kwon, H., & Kim, K. (2017). Approximate Subject Specific Pseudo MRI from an Available MRI Dataset for MEG Source Imaging. Frontiers in neuroinformatics, 11, 50.](https://www.frontiersin.org/articles/10.3389/fninf.2017.00050/full)
 
- [Whalen, C., Maclin, E. L., Fabiani, M., & Gratton, G. (2008). Validation of a method for coregistering scalp recording locations with 3D structural MR images. Human Brain Mapping, 29(11), 1288-1301.](http://onlinelibrary.wiley.com/doi/10.1002/hbm.20465/full)



