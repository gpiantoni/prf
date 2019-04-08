%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%
%
% Script to initialize the computation of pRF parameters from nifti files
% preprocessed with either UMCU or NYU pipline.
%
% Input: <subjectcode>
% Output: pRF parameters <ang>, <ecc>, <expt>, <rfsize>, <R2>, <gain> found
% in <output_dir>.
%


%%
% % % % 
% % % % % Read in niftis (preprocessed) 
% % % % img = 'sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii';
% % % % img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
% % % % run1 = niftiread (img);
% % % % img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
% % % % img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz';
% % % % run2 = niftiread (img);
% % % % 
% % % % % Speficy vector of voxel indices to analyze
% % % % vxs = zeros(size(run1, 1), size(run1, 2), size(run1, 3));
% % % % vxs(30:60, 30:60, 26) = run1(30:60, 30:60, 26) >= 0;
% % % 
% % % %%
% % % 
% % % % Compute pRF parameters with GLMDenoise
% % % results  = analyzePRF({images, images},{run1, run2}, TR, struct('seedmode',[0 1],'display','on', 'wantglmdenoise', 1, 'vxs', find(vxs)));

%%

addpath(genpath('/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF'))   %% TO REMOVE

addpath(genpath('/Fridge/users/margriet/projects/analysis/analyzeprf'))

parpool(40)


%%

% Enter subject code
subjectcode = 'sub-visual01';

%% NYU preprocessing pipeline

%%%%%%%%%%%% 3TMB %%%%%%%%%%%%


nifti =  {'/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz',
    '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz'};        % TO REMOVE



nifti_gio =  {['/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
    ['/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};    % TO REMOVE
    
nifti_margriet =  {['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
    ['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};


output_dir = '/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF/results_analyzePRF/01';     % TO REMOVE

output_dir_nyu3TMB = ['Fridge/users/margriet/projects/analysis/analyzeprf/results/nyu/', subjectcode, '/ses-UMCU3TMB'];


threshold = 1200;

compute_prf(nifti, output_dir_nyu3TMB, threshold)


%%%%%%%%%%%% 7TGE %%%%%%%%%%%%

nifti =  {'',
    ''};
output_dir_nyu7TGE = '';
threshold = 1200;

compute_prf(nifti, output_dir_nyu7TGE, threshold)





%%%%%%%%%%%% 7TSE %%%%%%%%%%%%



%% UMCU preprocessing pipeline

%%%%%%%%%%%% 3TMB %%%%%%%%%%%%

nifti =  {'',
    ''};
output_dir_umcu3TMB = '';
threshold = 1200;

compute_prf(nifti, output_dir_umcu3TMB, threshold)

%%%%%%%%%%%% 7TGE %%%%%%%%%%%%
nifti =  {'',
    ''};
output_dir_umcu7TGE = '';
threshold = 1200;

compute_prf(nifti, output_dir_umcu7TGE, threshold)



%%

nifti = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';

hdr = niftiinfo(nifti);


nii = niftiread(nifti);
mean_nii = double(mean(nii, 4));


hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';
niftiwrite(mean_nii, '/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF/results_analyzePRF/01/mean.nii', hdr);
