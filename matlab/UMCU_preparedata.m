% % =========  SCRIPT prepare UMCU data for analyzePRF  ========= % 
%
% Script (1) averages two pRF runs and (2) splits up merged pRF runs
%
% Input: <subjectcode>, <session>
% Output: <bairprf_AVERAGED_bold.nii>
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

%% Specify parameters

subjectcode = 'sub-visual06';              
session = 'ses-UMCU3TMB';
output_dir_averaged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session];
output_dir_merged = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold'];

%% Run functions

% ========= AVERAGE SEPARATE RUNS ========= %
average_runs (subjectcode, session, output_dir_averaged)
% ========= SPLIT MERGED RUNS ========= %
split_merged_runs (subjectcode, session, output_dir_merged)

disp ('Done')