function compute_prf(subjectcode, session, nifti, n_volumes, output_dir, threshold, AnalyzeMergedRuns, AnalyzeSeparateRuns, AnalyzeAveragedRuns, UseDenoise)
%
% function COMPUTE_PRF (nifti, output_dir, threshold)
%
% Inputs:
% <nifti> is the path to the nifti file, as string
% <n_volumes> is a list of the number of volumes for each run
% <output_dir> 
% <threshold> voxels with values greater than the given threshold will be
% included in the computation. 
% <AnalyzeMergedRuns> true/false
% <UseDenoise> true/false
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
%% 

disp('-- Loading stimuli --')

% ========= Read in parameters ========= % 

%% MERGED RUNS 
    % ========= MERGED RUNS ========= % 

if AnalyzeMergedRuns == true
   disp('-- Starting analyzePRF for merged runs --')
    hdr = niftiinfo(nifti);                           % read nifti header
    TR = hdr.PixelDimensions(4);                      % 850 ms
    if TR ~= 0.85
        % fprintf('TR from the nifti is %.fs but it should be 0.850s. Fixing it\n', TR)
        TR = 0.850;
    end    

    images = {};
    images{1} = read_bair_stimuli(subjectcode, session, n_volumes, TR);

    disp('-- Loading fMRI --')
    nii = {};
    nii{1} = niftiread(nifti);

    first_nii = nii{1}(:, :,:, 1);      % reference scan: first volume
    vxs = find(first_nii > threshold);  % find voxels above threshold

    n_dim = size(nii{1});
    n_vox = prod(n_dim(1:3));
    fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)
end

%% SEPARATE RUNS
    % ========= SEPARATE RUNS ========= % 
if AnalyzeSeparateRuns == true
    disp('-- Starting analyzePRF for separate runs --')
    % MAX_N_VOLUMES = 252;        % previously 248
    images = {};

    % 1st run
    run = 1;
    hdr = niftiinfo(nifti{1});                        % read nifti header
    TR = hdr.PixelDimensions(4);                      % 850 ms
    n_volumes = hdr.ImageSize(4);                     % no. of dynamics
    % n_volumes = min(n_volumes, MAX_N_VOLUMES);
    images_run1 = read_bair_stimuli(subjectcode, session, run, n_volumes, TR);

    % 2nd run
    run = 2;
    hdr = niftiinfo(nifti{2});                        % read nifti header
    TR = hdr.PixelDimensions(4);                      % 850 ms
    n_volumes = hdr.ImageSize(4);                     % no. of dynamics
    % n_volumes = min(n_volumes, MAX_N_VOLUMES);
    images_run2 = read_bair_stimuli(subjectcode, session, run, n_volumes, TR);
    
    images = {images_run1, images_run2};

% %     % Resize nifti to lowest no. of dynamics
% %     disp('-- Loading fMRI --')
% %     nii = {};
% %     for i = 1:length(nifti)                          % length(nifti) = 2
% %         temp = niftiread(nifti{i});
% %         n_volumes = size(temp, 4);
% %         n_volumes = min(n_volumes, MAX_N_VOLUMES);
% %         nii{i} = temp(:, :, :, 1:n_volumes);
% %     end

    nii = {};
    for i = 1:length(nifti)
        nii{i} = niftiread (nifti{i});
    end
       
    % Reference
    first_nii = nii{1}(:, :,:, 1);      % reference scan: first volume of first run
    vxs = find(first_nii > threshold);  % find voxels above threshold

    n_dim = size(nii{1});
    n_vox = prod(n_dim(1:3));
    fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)
end

%% AVERAGED RUNS
        % ========= AVERAGED RUNS ========= %
if AnalyzeAveragedRuns == true
   disp('-- Starting analyzePRF for averaged runs --')
    hdr = niftiinfo(nifti);                           % read nifti header
    TR = hdr.PixelDimensions(4);                      % 850 ms
    if TR ~= 0.85
        % fprintf('TR from the nifti is %.fs but it should be 0.850s. Fixing it\n', TR)
        TR = 0.850;
    end    
    
    images = {};
    images{1} = read_bair_stimuli(subjectcode, session, n_volumes, TR);

    disp('-- Loading fMRI --')
    nii = {};
    nii{1} = niftiread(nifti);

    first_nii = nii{1}(:, :,:, 1);      % reference scan: first volume
    vxs = find(first_nii > threshold);  % find voxels above threshold

    n_dim = size(nii{1});
    n_vox = prod(n_dim(1:3));
    fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)
end

%%   
 % ========= Compute pRF parameters WITH/WITHOUT GLMdenoise ========= % 

if UseDenoise == false
    disp('%%%%%%%%%%%% Running analyzePRF without GLMdenoise %%%%%%%%%%%%')
    results = analyzePRF(images, nii, TR, struct('seedmode',[0 1 2],'display','off', 'wantglmdenoise', 0, 'vxs', vxs));
else
    disp('%%%%%%%%%%%% Running analyzePRF with GLMdenoise %%%%%%%%%%%%')
    results = analyzePRF(images, nii, TR, struct('seedmode',[0 1 2],'display','off', 'wantglmdenoise', 1, 'vxs', vxs));
end

%%
% ========= OUTPUT ========= %

hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';

% Output parameters saved in <output_dir>
fields = {'ang', 'ecc', 'expt', 'rfsize', 'R2', 'gain'};
for i = 1:length(fields)
    niftiwrite(results.(fields{i}), fullfile(output_dir, [fields{i}  '.nii']), hdr);
end

%% End
