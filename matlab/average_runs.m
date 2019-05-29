% % =========  SCRIPT Compute average pRF timeseries  ========= % 
%
% Script averages the 2 pRF runs
%
% Input: <subjectcode>, <session>
% Output: <bairprf_AVERAGED_bold.nii>
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

%% Specify parameters

subjectcode = 'sub-visual06';              
subjectnumber = str2num(subjectcode (11:12));
session = 'ses-UMCU3TMB';

%%

% ========= Read in nifti  ========= % 
if session == 'ses-UMCU3TMB';
    merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_bold-rwm.nii']);
    run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-rwm.nii']);
    dyn_run1 = size(run1,4);
    run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-rwm.nii']);
    dyn_run2 = size(run2,4);    
    
    n_vol = min(dyn_run1, dyn_run2);
    
    part1 = merged (:,:,:, 1:n_vol);
    part2 = merged (:,:,:, (1+n_vol):(2*n_vol));
    
else
    merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_bold-masked-mc-warp.nii']);
    run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-masked-mc-warp.nii']);
    dyn_run1 = size(run1,4);
    run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold-masked-mc-warp.nii']);
    dyn_run2 = size(run2,4);
 
    n_vol = min(dyn_run1, dyn_run2);
    
    part1 = merged (:,:,:, 1:n_vol);
    part2 = merged (:,:,:, (1+n_vol):(2*n_vol));
    
end

% ========= Take average of 2 pRF runs  ========= % 
avrg = (part1+part2)/2;

figure(1)
plot (squeeze(part1(71,69,1, :)), 'LineWidth', 2, 'Color', 'r')
hold on 
plot (squeeze(part2(71,69,1, :)), 'LineWidth', 2, 'Color', 'm')
plot (squeeze(merged(71,69,1, :)), 'LineWidth', 2, 'Color', 'c')
plot (squeeze(avrg(71,69,1, :)), 'LineWidth', 2, 'Color', 'k')
legend ('Run1', 'Run2', 'Merged', 'Average')
% axis ([0 500 -20 20])
hold off

% ========= Write output nifti (task-bairprf_AVERAGED_bold)  ========= % 
path = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session];

if session == 'ses-UMCU3TMB'
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-rwm.nii']);
    filename = [subjectcode, '_', session, '_task-bairprf_AVERAGED_bold-rwm.nii'];
    niftiwrite (avrg, fullfile(path, filename), hdr)
else
    hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold-masked-mc-warp.nii']);
    filename = [subjectcode, '_', session, '_task-bairprf_AVERAGED_bold-masked-mc-warp.nii'];
    niftiwrite (avrg, fullfile(path, filename), hdr)
end
