%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%

%% Dataset: Visual01 3T

% Read in niftis (preprocessed) 
img = 'sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii';
img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
run1 = niftiread (img);
img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz';
run2 = niftiread (img);

% Speficy vector of voxel indices to analyze
vxs = zeros(size(run1, 1), size(run1, 2), size(run1, 3));
vxs(30:60, 30:60, 26) = run1(30:60, 30:60, 26) >= 0;

%%

% Compute pRF parameters with GLMDenoise
results  = analyzePRF({images, images},{run1, run2}, TR, struct('seedmode',[0 1],'display','on', 'wantglmdenoise', 1, 'vxs', find(vxs)));

%%

addpath(genpath('/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF'))



%%%%%% Visual01
%%%%  3TMB
% NYU preproc
nifti =  {'/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz',
    '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz'};
output_dir = '/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF/results_analyzePRF/01';
threshold = 1200;

compute_prf(nifti, output_dir, threshold)

% Wouter preproc
nifti =  {'',
    ''};
output_dir = '';
threshold = 1200;

compute_prf(nifti, output_dir, threshold)

%%%%  7TGE
nifti =  {'',
    ''};
output_dir = '';
threshold = 1200;

compute_prf(nifti, output_dir, threshold)
