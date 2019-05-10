% =========  SCRIPT analyzePRF UMCU data  ========= % 
%
% Script to initialize the computation of pRF parameters from nifti files
% preprocessed with either UMCU or NYU pipeline.
%
% Input: <subjectcode>
% Output: pRF parameters <ang>, <ecc>, <expt>, <rfsize>, <R2>, <gain>, <meanvol> found
% in <output_dir>.
%
%%
clear;
addpath(genpath('/Fridge/users/margriet/projects/analysis/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))

%% Specify subject code

subjectcode = 'sub-visual04';               % Enter subject code
subjectnumber = str2num(subjectcode (11:12));

Analyze3TMB = 0;
Analyze7TGE = 1;
Analyze7TSE = 0;

%% Log thresholds
 
% =========  Threshold UMCU preprocessed data ========= % 
thresholdUMCU_3TMB = [700, 750, 200, 550, 800, 950  NaN, 600, NaN, NaN, 450, 250];
thresholdUMCU_7TGE = [150, 100, 100, 100, 100, 300, NaN, 250, 250, 250, 200, 150];
thresholdUMCU_7TSE = [50, 50, 50, NaN, 175, NaN, NaN, 50, NaN, 100, 300];

% =========  Threshold NYU preprocessed data ========= % 
thresholdNYU_3TMB = [800, 700];


%% Start parallel pool
% parpool(40)       % Parallel pool

%%  Merged pRF runs

% =========  3T (MB) ========= % 

if Analyze3TMB == true
    disp('%%%%%%%%%%%% MERGED NIFTIS - Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')

    nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold.nii']);
    nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold.nii']);
    n_volumes_run1 = size(nifti_run1, 4);
    n_volumes_run2 = size(nifti_run2, 4);
    n_volumes = [n_volumes_run1, n_volumes_run2];

    nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold-rwm.nii'];

    output_dir_umcu3TMB_merged = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/merged_bairprf'];

    threshold = thresholdUMCU_3TMB(subjectnumber);         

    compute_prf(nifti_merged, n_volumes, output_dir_umcu3TMB_merged, threshold)
end

% =========  7T (GE) ========= % 

if Analyze7TGE == true
    disp('%%%%%%%%%%%% MERGED NIFTIS - Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')

    nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold.nii']);
    nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold.nii']);
    n_volumes_run1 = size(nifti_run1, 4);
    n_volumes_run2 = size(nifti_run2, 4);
    n_volumes = [n_volumes_run1, n_volumes_run2];

    nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold-masked-mc-warp.nii'];

    output_dir_umcu7TGE_merged = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/merged_bairprf'];

    threshold = thresholdUMCU_7TGE(subjectnumber);        

    compute_prf(nifti_merged, n_volumes, output_dir_umcu7TGE_merged, threshold)
end

% =========  7T (SE) ========= % 

if Analyze7TSE == true
    disp('%%%%%%%%%%%% MERGED NIFTIS -  Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

    nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold.nii']);
    nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold.nii']);
    n_volumes_run1 = size(nifti_run1, 4);
    n_volumes_run2 = size(nifti_run2, 4);
    n_volumes = [n_volumes_run1, n_volumes_run2];

    nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold-masked-mc-warp.nii'];

    output_dir_umcu7TSE_merged = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/merged_bairprf'];

    threshold = thresholdUMCU_7TSE(subjectnumber);        

    compute_prf(nifti_merged, n_volumes, output_dir_umcu7TSE_merged, threshold)
end


%% UMCU preprocessing pipeline (3T, 7TGE, 7TSE) - Analyze pRF runs separately

% % % =========  3T (MB) ========= % 
% % 
% % if Analyze3TMB == true
% %     disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')
% % 
% %     nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii'],
% %         ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-02_bold-rwm.nii']};
% % 
% %     output_dir_umcu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB'];
% % 
% %     % % % % % Run with GLM denoise
% %     % % output_dir_umcu3TMB_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU3TMB'];
% %     % % output_dir_umcu3TMB = output_dir_umcu3TMB_denoise;
% % 
% %     threshold = thresholdUMCU_3TMB(subjectnumber);         
% % 
% %     compute_prf(nifti, output_dir_umcu3TMB, threshold)
% % end
% % 
% % % =========  7T (GE) ========= % 
% % 
% % if Analyze7TGE == true
% %     disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')
% % 
% %     nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-01_bold-masked-mc-warp.nii'],
% %         ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-02_bold-masked-mc-warp.nii']};
% % 
% %     output_dir_umcu7TGE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE'];
% % 
% % 
% %     % % % % % % Run with GLM denoise
% %     % % output_dir_umcu7TGE_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TGE'];
% %     % % output_dir_umcu7TGE = output_dir_umcu7TGE_denoise;
% % 
% %     threshold = thresholdUMCU_7TGE(subjectnumber);         
% % 
% %     compute_prf(nifti, output_dir_umcu7TGE, threshold)
% % end
% %      
% % % =========  7T (SE) ========= %     
% % 
% % if Analyze7TSE == true
% %     disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')
% % 
% %     nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-01_bold-masked-mc-warp.nii'],
% %         ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-02_bold-masked-mc-warp.nii']};
% % 
% %     output_dir_umcu7TSE = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE'];
% % 
% %     % % % % % % Run with GLM denoise
% %     % % output_dir_umcu7TSE_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TSE'];
% %     % % output_dir_umcu7TSE = output_dir_umcu7TSE_denoise;
% % 
% %     threshold = thresholdUMCU_7TSE(subjectnumber);       
% % 
% %     compute_prf(nifti, output_dir_umcu7TSE, threshold)
% % end


%% NYU preprocessing pipeline (3T)
% % 
% % % =========  3T (MB) ========= % 
% % 
% % if Analyze3TMB == true
% %     disp('%%%%%%%%%%%% Starting analyzePRF for NYU 3T (MB) %%%%%%%%%%%%')
% % 
% %     nifti =  {['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
% %         ['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};
% % 
% %     output_dir_nyu3TMB = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results/nyu/', subjectcode, '/ses-UMCU3TMB'];
% %     output_dir_nyu3TMB_denoise = ['/Fridge/users/margriet/projects/analysis/analyzeprf/results_glmdenoise/nyu/', subjectcode, '/ses-UMCU3TMB'];
% % 
% %     % % % Run with GLM denoise
% %     % % output_dir_nyu3TMB = output_dir_nyu3TMB_denoise;
% % 
% %     threshold = thresholdNYU_3TMB(subjectnumber);       
% % 
% %     compute_prf(nifti, output_dir_nyu3TMB, threshold)
% % end

%% Manually create R2 mask 

% % % fslmaths R2.nii -nan -thr 5 -bin R2mask_5.nii
% % % fslmaths R2.nii -nan -thr 10 -bin R2mask_10.nii
% % % fslmaths R2.nii -nan -thr 15 -bin R2mask_15.nii


