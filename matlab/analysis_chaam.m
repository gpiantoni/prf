% Script analyzePRF sub-chaam

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

%% Specify parameters

subjectcode = 'sub-chaam';              
session = 'ses-UMCU3Tday139';
threshold3T = 500;
threshold7T = 75;

UseDenoise = false;

Analyze3TMB = 1;
Analyze7TGE = 0;
Analyze7TSE = 0;

AnalyzeMergedRuns = false;
AnalyzeSeparateRuns = true;
AnalyzeAveragedRuns = false;

%% Split up merged runs

if session == 'ses-UMCU3Tday139'
    merged = niftiread ('/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_bold-rwm.nii');
    
    run1 = niftiread ('/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-1_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-1_bold-rwm.nii');
    dyn_run1 = size(run1,4);
    run2 = niftiread ('/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-2_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-2_bold-rwm.nii');
    dyn_run2 = size(run2,4);
    n_vol = min(dyn_run1, dyn_run2);
    
    merged_part1 = merged (:,:,:, 1:dyn_run1);
    merged_part2 = merged (:,:,:, (1+dyn_run1):(dyn_run1+dyn_run2));
end
   
% ========= Write output nifti (task-bairprf_MERGED_runxx_bold)  ========= % 
output_dir = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold'];

if session == 'ses-UMCU3Tday139'
    % Run1
    hdr = niftiinfo ('/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-1_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-1_bold-rwm.nii');
    filename_run1 = [subjectcode, '_', session, '_task-bairprf_MERGED_run01_bold-rwm.nii'];
    niftiwrite (merged_part1, fullfile(output_dir, filename_run1), hdr)
    hdr = niftiinfo ('/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-2_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_run-2_bold-rwm.nii');
    filename_run2 = [subjectcode, '_', session, '_task-bairprf_MERGED_run02_bold-rwm.nii'];
    niftiwrite (merged_part2, fullfile(output_dir, filename_run2), hdr)
end
   
disp (['Done splitting merged pRF runs for ', subjectcode, ': ', session])

%%

if AnalyzeSeparateRuns == true
         disp('-- Starting analyzePRF for separate runs - new model --')
    % =========  3T (MB) ========= % 

    if Analyze3TMB == true
        session =  'ses-UMCU3Tday139';
        disp('%%%%%%%%%%%% Starting analyzePRF for UMCU 3T (MB) %%%%%%%%%%%%')
       
        % Read in split up niftis
        nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_run01_bold-rwm.nii'],
             ['/Fridge/users/margriet/subjects/bids_umcupreproc/sub-chaam/ses-UMCU3Tday139/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_bold/sub-chaam_ses-UMCU3Tday139_task-bairprf_MERGED_run02_bold-rwm.nii']};
          
        % Read no. of dynamics per run
        nifti_run1 = niftiread (nifti{1});
        nifti_run2 = niftiread (nifti{2});
        n_volumes_run1 = size(nifti_run1, 4);
        n_volumes_run2 = size(nifti_run2, 4);
        n_volumes = [n_volumes_run1, n_volumes_run2];

        threshold3T = 500;    

        output_dir_umcu3T = ['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/sub-chaam/ses-UMCU3Tday139/separate_bairprf'];
        
        % ========= ANALYZEPRF ========= % 
        compute_prf_patient(subjectcode, session, nifti, n_volumes, output_dir_umcu3T, threshold3T, AnalyzeMergedRuns, AnalyzeSeparateRuns, AnalyzeAveragedRuns, UseDenoise)
        % ======= PIXELS2DEGREES ======= % 
        convert_pixels2degrees_patient (session, output_dir_umcu3T)
        % ======== MEAN VOLUME ========= % 
        compute_mean_volume_patient (subjectcode, session, output_dir_umcu3T)
    end
end

                

