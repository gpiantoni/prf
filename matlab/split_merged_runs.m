% % =========  SCRIPT split up merged pRF timeseries  ========= % 
%
% Script splits the merged nifti up in 2 separate pRF runs
%
% Input: <subjectcode>, <session>
% Output: <bairprf_AVERAGED_bold.nii>
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

%% Specify parameters

subjectcode = 'sub-visual03';              
session = 'ses-UMCU7TSE';

%%

% ========= Read in nifti  ========= % 
if session == 'ses-UMCU3TMB'
    merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_bold-rwm.nii']);
    run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-rwm.nii']);
    dyn_run1 = size(run1,4);
    run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-rwm.nii']);
    dyn_run2 = size(run2,4);    
    
    n_vol = min(dyn_run1, dyn_run2);
    
    merged_part1 = merged (:,:,:, 1:dyn_run1);
    merged_part2 = merged (:,:,:, (1+dyn_run1):(dyn_run1+dyn_run2));
    
else
    merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_bold-masked-mc-warp.nii']);
    run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-masked-mc-warp.nii']);
    dyn_run1 = size(run1,4);
    run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-masked-mc-warp.nii']);
    dyn_run2 = size(run2,4);
 
    n_vol = min(dyn_run1, dyn_run2);
    
    merged_part1 = merged (:,:,:, 1:dyn_run1);
    merged_part2 = merged (:,:,:, (1+dyn_run1):(dyn_run1+dyn_run2));
end

% ========= Write output nifti (task-bairprf_MERGED_runxx_bold)  ========= % 
path = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold'];

if session == 'ses-UMCU3TMB'
    % Run1
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-rwm.nii']);
    filename_run1 = [subjectcode, '_', session, '_task-bairprf_MERGED_run01_bold-rwm.nii'];
    niftiwrite (merged_part1, fullfile(path, filename_run1), hdr)
    % Run2
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-rwm.nii']);
    filename_run2 = [subjectcode, '_', session, '_task-bairprf_MERGED_run02_bold-rwm.nii'];
    niftiwrite (merged_part2, fullfile(path, filename_run2), hdr)
   
else
    % Run1
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-masked-mc-warp.nii']);
    filename_run1 = [subjectcode, '_', session, '_task-bairprf_MERGED_run01_bold-masked-mc-warp.nii'];
    niftiwrite (merged_part1, fullfile(path, filename_run1), hdr)
    % Run2
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-masked-mc-warp.nii']);
    filename_run2 = [subjectcode, '_', session, '_task-bairprf_MERGED_run02_bold-masked-mc-warp.nii'];
    niftiwrite (merged_part2, fullfile(path, filename_run2), hdr)
end

%%

% figure(1)
% plot (squeeze(merged_part1(71,69,1, :)), 'LineWidth', 2, 'Color', 'r')
% hold on 
% plot (squeeze(merged_part2(71,69,1, :)), 'LineWidth', 2, 'Color', 'm')
% plot (squeeze(merged(71,69,1, :)), 'LineWidth', 2, 'Color', 'c')
% legend ('Run1', 'Run2', 'Merged')
% % axis ([0 500 -20 20])
% hold off


disp ('Done')















