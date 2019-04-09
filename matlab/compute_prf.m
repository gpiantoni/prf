function compute_prf(nifti, output_dir, threshold)
%
% function COMPUTE_PRF (nifti, output_dir, threshold)
%
% Inputs:
% <nifti> provides the data as a cell vector of voxels x time.  can also be
% X x Y x Z x time.
% <output_dir> 
% <threshold> voxels with values greater than the given threshold will be
% included in the computation. 
%
% Outputs:
% <ang> contains pRF angle estimates.  Values range between 0 and 360 degrees.
%   0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical
%   meridian, and so on.
% <ecc> contains pRF eccentricity estimates.  Values are in pixel units with a lower
%   bound of 0 pixels.
% <rfsize> contains pRF size estimates.  pRF size is defined as sigma/sqrt(n) where
%   sigma is the standard of the 2D Gaussian and n is the exponent of the power-law
%   function.  Values are in pixel units with a lower bound of 0 pixels.
% <expt> contains pRF exponent estimates.
% <gain> contains pRF gain estimates.  Values are in the same units of the data
%   and are constrained to be non-negative.
% <R2> contains R^2 values that indicate the goodness-of-fit of the model to the data.
%   Values are in percentages and generally range between 0% and 100%.  The R^2 values
%   are computed after projecting out polynomials from both the data and the model fit.
%   (Because of this projection, R^2 values can sometimes drop below 0%.)

hdr = niftiinfo(nifti{1});          % read nifti header
n_volumes = hdr.ImageSize(4);       % no. of dynamics
TR = hdr.PixelDimensions(4);        % 850 ms

disp('Loading stimuli')
images = read_bair_stimuli(n_volumes, TR);

disp('Loading images')
for i = 1:length(nifti)
    nii{i} = niftiread(nifti{i});
end

first_nii = nii{1}(:, :,:, 1);      % reference scan: first volume
vxs = find(first_nii > threshold);  % find voxels above threshold

n_dim = size(nii{1});
n_vox = prod(n_dim(1:3));
fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)

% Compute pRF parameters without GLMdenoise
results = analyzePRF({images, images}, nii, TR, struct('seedmode',[0 1 2],'display','off', 'wantglmdenoise', 0, 'vxs', vxs));

hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';

% Output parameters saved in <output_dir>
fields = {'ang', 'ecc', 'expt', 'rfsize', 'R2', 'gain'};
for i = 1:length(fields)
    niftiwrite(results.(fields{i}), fullfile(output_dir, [fields{i}  '.nii']), hdr);
end

