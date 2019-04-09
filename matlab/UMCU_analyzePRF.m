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

subjectcode = 'sub-visual02';               % Enter subject code
number = str2num(subjectcode (11:12));

%% Log thresholds

% Threshold UMCU preprocessed data
thresholdUMCU_3TMB = [500, 750];
thresholdUMCU_7TGE = [150, 100];
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
threshold = thresholdUMCU_3TMB(number);         

compute_prf(nifti, output_dir_umcu3TMB, threshold)


%%%%%%%%%%%% 7T (GE) %%%%%%%%%%%%

disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-01_bold-rwm.nii'],
    ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-02_bold-rwm.nii']};

output_dir_umcu7TGE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE'];
threshold = thresholdUMCU_7TGE(number);       

compute_prf(nifti, output_dir_umcu7TGE, threshold)
    
%%%%%%%%%%%% 7T (SE) %%%%%%%%%%%%        % Not preprocessed yet
%
% disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

%% NYU preprocessing pipeline

%%%%%%%%%%%% 3T (MB) %%%%%%%%%%%%

disp('%%%%%%%%%%%% Starting analyzePRF for NYU 3T (MB) %%%%%%%%%%%%')

nifti =  {['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
    ['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};

output_dir_nyu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/nyu/', subjectcode, '/ses-UMCU3TMB'];
threshold = thresholdNYU_3TMB(number);        

compute_prf(nifti, output_dir_nyu3TMB, threshold)

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
%
% % nifti = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
% % 
% % hdr = niftiinfo(nifti);
% % 
% % nii = niftiread(nifti);
% % mean_nii = double(mean(nii, 4));
% % 
% % hdr.ImageSize = hdr.ImageSize(1:3);
% % hdr.PixelDimensions = hdr.PixelDimensions(1:3);
% % hdr.Datatype = 'double';
% % niftiwrite(mean_nii, '/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF/results_analyzePRF/01/mean.nii', hdr);
% % 
% %  



