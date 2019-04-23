%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%
%
% Script to initialize the computation of pRF parameters from nifti files
% preprocessed with either UMCU or NYU pipeline.
%
% Input: <subjectcode>
% Output: pRF parameters <ang>, <ecc>, <expt>, <rfsize>, <R2>, <gain> found
% in <output_dir>.
%

%% Specify subject code

subjectcode = 'sub-visual03';               % Enter subject code
number = str2num(subjectcode (11:12));

%% Log thresholds

% Threshold UMCU preprocessed data
thresholdUMCU_3TMB = [500, 750, 800, NaN, 700, 950  NaN, NaN, NaN];
thresholdUMCU_7TGE = [150, 100, 100, NaN, 200, 300, NaN, NaN, 250];
thresholdUMCU_7TSE = [];

% Threshold NYU preprocessed data
thresholdNYU_3TMB = [800, 700];
thresholdNYU_7TGE = [];
thresholdNYU_7TSE = [];

%% Start parallel pool

addpath(genpath('/Fridge/users/margriet/projects/analysis/analyzeprf'))
% parpool(40)       % Parallel pool


%% UMCU preprocessing pipeline

%%%%%%%%%%%% 3T (MB) %%%%%%%%%%%%

disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii'],
    ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-02_bold-rwm.nii']};

output_dir_umcu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB'];
% output_dir_umcu3TMB_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU3TMB'];

% % % % % Run with GLM denoise
% % % % output_dir_umcu3TMB = output_dir_umcu3TMB_denoise;

threshold = thresholdUMCU_3TMB(number);         

compute_prf(nifti, output_dir_umcu3TMB, threshold)


%%%%%%%%%%%% 7T (GE) %%%%%%%%%%%%
% % % 
% % % disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')
% % % 
% % % nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-01_bold-rwm.nii'],
% % %     ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-02_bold-rwm.nii']};
% % % 
% % % output_dir_umcu7TGE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE'];
% % % % output_dir_umcu7TGE_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TGE'];
% % % 
% % % % % % % % % Run with GLM denoise
% % % % % % % % output_dir_umcu7TGE = output_dir_umcu7TGE_denoise;
% % % 
% % % threshold = thresholdUMCU_7TGE(number);       
% % % 
% % % compute_prf(nifti, output_dir_umcu7TGE, threshold)

    
%%%%%%%%%%%% 7T (SE) %%%%%%%%%%%%        % Not preprocessed yet
%
% disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

%% NYU preprocessing pipeline

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
% % % % threshold = thresholdNYU_3TMB(number);       
% % % % 
% % % % compute_prf(nifti, output_dir_nyu3TMB, threshold)

% %%%%%%%%%%%% 7T (GE) %%%%%%%%%%%%      % Not preprocessed yet
%
% disp('%%%%%%%%%%%% Starting analyzePRF for NYU 7T (GE) %%%%%%%%%%%%')
%
% nifti =  {'',
%     ''};
% output_dir_nyu7TGE = '';
% Enter threshold 
% threshold = thresholdNYU_7TGE(number);
% 
% compute_prf(nifti, output_dir_nyu7TGE, threshold)


%%%%%%%%%%%% 7T (SE) %%%%%%%%%%%%        % Not preprocessed yet
%
% disp('%%%%%%%%%%%% Starting analyzePRF for NYU 7T (SE) %%%%%%%%%%%%')


%% Manually compute mean

% % % Enter preprocessed nifti
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
% % 
% % 
% % 
% % 








