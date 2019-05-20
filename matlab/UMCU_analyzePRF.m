% % =========  SCRIPT analyzePRF UMCU data  ========= % 
%
% Script to initialize the computation of pRF parameters from nifti files
% preprocessed with either UMCU or NYU pipeline.
%
% Input: <subjectcode>, <UseDenoise>, <Analyze3TMB>, <Analyze7TGE>,
% <Analyze7TSE>, <AnalyzeMergedRuns>
% Output: pRF parameters <ang>, <ecc>, <ecc_deg>, <expt>, <rfsize>, <rfsize_deg>
% <R2>, <gain>, <mean> found in <output_dir>.

%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))

%% Specify parameters

subjectcode = 'sub-visual06';               
subjectnumber = str2num(subjectcode (11:12));

UseDenoise = false;

Analyze3TMB = 0;
Analyze7TGE = 1;
Analyze7TSE = 0;

AnalyzeMergedRuns = false;

%% Log thresholds
 
% =========  Threshold UMCU preprocessed data ========= % 
thresholdUMCU_3TMB = [700, 750, 200, 550, 800, 950  NaN, 600, NaN, NaN, 450, 250];
thresholdUMCU_7TGE = [150, 100, 100, 100, 100, 300, NaN, 250, 250, 250, 200, 150];
thresholdUMCU_7TSE = [50, 50, 50, NaN, 175, 50, NaN, 50, NaN, 100, 300, 120];

% =========  Threshold NYU preprocessed data ========= % 
thresholdNYU_3TMB = [800, 700];


%% Start parallel pool
% parpool(40)       % Parallel pool

%%  Merged pRF runs

if AnalyzeMergedRuns == true

    % =========  3T (MB) ========= % 

    if Analyze3TMB == true
        session = 'ses-UMCU3TMB';
        disp('%%%%%%%%%%%% MERGED NIFTIS - Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')

        nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold.nii']);
        nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold.nii']);
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];

        nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold-rwm.nii'];

        output_dir_umcu3TMB_merged = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/merged_bairprf'];
        output_dir_umcu3TMB_merged_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU3TMB/merged_bairprf'];

        threshold = thresholdUMCU_3TMB(subjectnumber);

            if UseDenoise == true
                output_dir_umcu3TMB_merged = output_dir_umcu3TMB_merged_denoise;
            end

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti_merged, n_volumes, output_dir_umcu3TMB_merged, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= % 
        convert_pixels2degrees (subjectcode, session, output_dir_umcu3TMB_merged)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu3TMB_merged)
    end

    % =========  7T (GE) ========= % 

    if Analyze7TGE == true
        session = 'ses-UMCU7TGE';
        disp('%%%%%%%%%%%% MERGED NIFTIS - Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')

        nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold.nii']);
        nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold.nii']);
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];

        nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold-masked-mc-warp.nii'];

        output_dir_umcu7TGE_merged = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/merged_bairprf'];
        output_dir_umcu7TGE_merged_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TGE/merged_bairprf'];

        threshold = thresholdUMCU_7TGE(subjectnumber);  

            if UseDenoise == true
                output_dir_umcu7TGE_merged = output_dir_umcu7TGE_merged_denoise;
            end

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti_merged, n_volumes, output_dir_umcu7TGE_merged, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= % 
        convert_pixels2degrees (subjectcode, session, output_dir_umcu7TGE_merged)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu7TGE_merged)
    end

    % =========  7T (SE) ========= % 

    if Analyze7TSE == true
        session = 'ses-UMCU7TSE';
        disp('%%%%%%%%%%%% MERGED NIFTIS -  Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

        nifti_run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold.nii']);
        nifti_run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold.nii']);
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];

        nifti_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold-masked-mc-warp.nii'];
        output_dir_umcu7TSE_merged = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/merged_bairprf'];
        output_dir_umcu7TSE_merged_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TSE/merged_bairprf'];

        threshold = thresholdUMCU_7TSE(subjectnumber); 

            if UseDenoise == true
                output_dir_umcu7TSE_merged = output_dir_umcu7TSE_merged_denoise;
            end

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti_merged, n_volumes, output_dir_umcu7TSE_merged, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= %  
        convert_pixels2degrees (subjectcode, session, output_dir_umcu7TSE_merged)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu7TSE_merged)       
    end
end

%% UMCU preprocessing pipeline (3T, 7TGE, 7TSE) - Analyze pRF runs separately

if AnalyzeMergedRuns == false
         disp('Starting analyzePRF for separate runs')
    % =========  3T (MB) ========= % 

    if Analyze3TMB == true
        session = 'ses-UMCU3TMB';
        disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')

        nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU3TMB_task-bairprf_run-02_bold-rwm.nii']};
        
        nifti_run1 = niftiread (nifti{1});
        nifti_run2 = niftiread (nifti{2});
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];

        threshold = thresholdUMCU_3TMB(subjectnumber);    

        output_dir_umcu3TMB = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB'];
        output_dir_umcu3TMB_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU3TMB'];

        if UseDenoise == true
            output_dir_umcu3TMB = output_dir_umcu3TMB_denoise;
        end

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti, n_volumes, output_dir_umcu3TMB, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= %  
        convert_pixels2degrees (subjectcode, session, output_dir_umcu3TMB)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu3TMB)  
    end

    % =========  7T (GE) ========= % 

    if Analyze7TGE == true
        session = 'ses-UMCU7TGE';
        disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (GE) %%%%%%%%%%%%')

        nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-01_bold-masked-mc-warp.nii'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TGE_task-bairprf_run-02_bold-masked-mc-warp.nii']};
        
        nifti_run1 = niftiread (nifti{1});
        nifti_run2 = niftiread (nifti{2});
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];
        
        threshold = thresholdUMCU_7TGE(subjectnumber);  

        output_dir_umcu7TGE = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE'];
        output_dir_umcu7TGE_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TGE'];

        if UseDenoise == true
            output_dir_umcu7TGE = output_dir_umcu7TGE_denoise;
        end   

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti, n_volumes, output_dir_umcu7TGE, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= %  
        convert_pixels2degrees (subjectcode, session, output_dir_umcu7TGE)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu7TGE)  
    end

    % =========  7T (SE) ========= %     

    if Analyze7TSE == true
        session = 'ses-UMCU7TSE';
        disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 7T (SE) %%%%%%%%%%%%')

        nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-01_bold-masked-mc-warp.nii'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02_bold/', subjectcode '_ses-UMCU7TSE_task-bairprf_run-02_bold-masked-mc-warp.nii']};
        
        nifti_run1 = niftiread (nifti{1});
        nifti_run2 = niftiread (nifti{2});
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];
        
        threshold = thresholdUMCU_7TSE(subjectnumber); 

        output_dir_umcu7TSE = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE'];
        output_dir_umcu7TSE_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/ses-UMCU7TSE'];

        if UseDenoise == true
            output_dir_umcu7TSE = output_dir_umcu7TSE_denoise;
        end  

        % ========= ANALYZEPRF ========= % 
        compute_prf(nifti, n_volumes, output_dir_umcu7TSE, threshold, AnalyzeMergedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= %  
        convert_pixels2degrees (subjectcode, session, output_dir_umcu7TSE)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume (subjectcode, session, output_dir_umcu7TSE)    
    end
end

%% NYU preprocessing pipeline (3T)
% % 
% % % =========  3T (MB) ========= % 
% % 
% % if Analyze3TMB == true
% %     session = 'ses-UMCU3TMB';
% %     disp('%%%%%%%%%%%% Starting analyzePRF for NYU 3T (MB) %%%%%%%%%%%%')
% % 
% %     nifti =  {['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
% %         ['/Fridge/users/margriet/subjects/bids_nyupreproc/data/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/volume/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']};
% % 
% %     threshold = thresholdNYU_3TMB(subjectnumber);  
% %     
% %     output_dir_nyu3TMB = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/nyu/', subjectcode, '/ses-UMCU3TMB'];
% %     output_dir_nyu3TMB_denoise = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/nyu/', subjectcode, '/ses-UMCU3TMB'];
% %    
% %     if UseDenoise == true
% %         output_dir_nyu3TMB = output_dir_nyu3TMB_denoise;
% %     end
% %     
% %     % ========= ANALYZEPRF ========= % 
% %     compute_prf(nifti_merged, n_volumes, output_dir_nyu3TMB, threshold, UseDenoise)
% %     % ======= PIXELS2DEGREES ======= %  
% %     convert_pixels2degrees (subjectcode, session, output_dir_nyu3TMB)
% %     % ======== MEAN VOLUME ========= % 
% %     compute_mean_volume (subjectcode, session, output_dir_nyu3TMB)  
% % end

%% Manually create R2 mask 

% % % fslmaths R2.nii -nan -thr 5 -bin R2mask_5.nii
% % % fslmaths R2.nii -nan -thr 10 -bin R2mask_10.nii
% % % fslmaths R2.nii -nan -thr 15 -bin R2mask_15.nii


%% End
disp ('end')

