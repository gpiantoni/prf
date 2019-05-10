% =========  SCRIPT Compute mean volume  ========= % 
%
% Script compute the mean volume of the merged time series.
%
% Input: <subjectcode>, <session>
% Output: <mean.nii> found in <output_dir>.
%

%%
addpath(genpath('/Fridge/users/margriet/projects/analysis/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))

%% Specify subject code

subject = 'sub-visual12';               % Enter subject code
session = 'ses-UMCU3TMB';

%% Compute mean

if session == 'ses-UMCU3TMB'
    nifti = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subject, '/', session, '/', subject, '_', session, '_task-bairprf_MERGED_bold/', subject, '_', session, '_task-bairprf_MERGED_bold-rwm'];
else 
    nifti = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subject, '/', session, '/', subject, '_', session, '_task-bairprf_MERGED_bold/', subject, '_', session, '_task-bairprf_MERGED_bold-masked-mc-warp.nii'];
end

hdr = niftiinfo(nifti);

nii = niftiread(nifti);
mean_nii = double(mean(nii, 4));

hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';

% ========= Write new nifti file ========= % 
output_dir = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subject, '/', session, '/merged_bairprf/mean.nii'];
niftiwrite(mean_nii, output_dir, hdr);

%%
