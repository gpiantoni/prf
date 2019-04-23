%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%
%
% Script to initialize the computation of pRF parameters from nifti files
% preprocessed with either UMCU or NYU pipeline.
%
% Input: <subjectcode>
% Output: pRF parameters <ang>, <ecc>, <expt>, <rfsize>, <R2>, <gain> found
% in <output_dir>.
%
%%
clear;
addpath(genpath('/Fridge/users/margriet/projects/analysis/analyzeprf'))


%% Specify subject code

subjectcode = 'sub-visual01';               % Enter subject code
subjectnumber = str2num(subjectcode (11:12));

%% Log thresholds

% Threshold UMCU preprocessed data
thresholdUMCU_3TMB = [500, 750, 800, NaN, 700, 950  NaN, NaN, NaN];
thresholdUMCU_7TGE = [150, 100, 100, NaN, 200, 300, NaN, NaN, 250];
thresholdUMCU_7TSE = [];

% Threshold NYU preprocessed data
thresholdNYU_3TMB = [800, 700];


%% Start parallel pool

% parpool(40)       % Parallel pool

%% UMCU preprocessing pipeline (3T, 7TGE, 7TSE)


%%%%%%%%%%%% 3T (MB) %%%%%%%%%%%%

disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii'],
    ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-02_bold-rwm.nii']};

output_dir_umcu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB'];

% % % % % Run with GLM denoise
% % output_dir_umcu3TMB_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU3TMB'];
% % output_dir_umcu3TMB = output_dir_umcu3TMB_denoise;

threshold = thresholdUMCU_3TMB(subjectnumber);         

compute_prf(nifti, output_dir_umcu3TMB, threshold)


%%%%%%%%%%%% 7T (GE) %%%%%%%%%%%%

disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-01_bold-rwm.nii'],
    ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-02_bold-rwm.nii']};

output_dir_umcu7TGE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE'];


% % % % % % Run with GLM denoise
% % output_dir_umcu7TGE_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TGE'];
% % output_dir_umcu7TGE = output_dir_umcu7TGE_denoise;

threshold = thresholdUMCU_7TGE(subjectnumber);         

compute_prf(nifti, output_dir_umcu7TGE, threshold)

    
%%%%%%%%%%%% 7T (SE) %%%%%%%%%%%%    

disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-01_bold-rwm.nii'],
    ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-02_bold-rwm.nii']};

output_dir_umcu7TSE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE'];

% % % % % % Run with GLM denoise
% % output_dir_umcu7TSE_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TSE'];
% % output_dir_umcu7TSE = output_dir_umcu7TSE_denoise;

threshold = thresholdUMCU_7TSE(subjectnumber);       

compute_prf(nifti, output_dir_umcu7TSE, threshold)


%% NYU preprocessing pipeline (3T)

%%%%%%%%%%%% 3T (MB) %%%%%%%%%%%%
% % % % 
% % % % disp('%%%%%%%%%%%% Starting analyzePRF for NYU 3T (MB) %%%%%%%%%%%%')
% % % % 
% % % % nifti =  {['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
% % % %     ['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};
% % % % 
% % % % output_dir_nyu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/nyu/', subjectcode, '/ses-UMCU3TMB'];
% % % % output_dir_nyu3TMB_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/nyu/', subjectcode, '/ses-UMCU3TMB'];
% % % % 
% % % % % % % Run with GLM denoise
% % % % % % output_dir_nyu3TMB = output_dir_nyu3TMB_denoise;
% % % % 
% % % % threshold = thresholdNYU_3TMB(subjectnumber);       
% % % % 
% % % % compute_prf(nifti, output_dir_nyu3TMB, threshold)


%% Manually compute mean

% % % Enter preprocessed nifti      (-rwm.nii)
% % nifti = '/Fridge/users/margriet/subjects/bids_umcupreproc/sub-visual03/ses-UMCU7TGE/sub-visual03_ses-UMCU7TGE_task-bairprf_run-01_bold/sub-visual03_ses-UMCU7TGE_task-bairprf_run-01_bold-rwm.nii';
% % 
% % hdr = niftiinfo(nifti);
% % 
% % nii = niftiread(nifti);
% % mean_nii = double(mean(nii, 4));
% % 
% % hdr.ImageSize = hdr.ImageSize(1:3);
% % hdr.PixelDimensions = hdr.PixelDimensions(1:3);
% % hdr.Datatype = 'double';
% % 
% % % Specify output_dir/mean.nii
% % niftiwrite(mean_nii, '/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/sub-visual03/ses-UMCU7TGE/mean.nii', hdr);









