function average_runs2 (subjectcode, session, output_dir)
%
% function AVERAGE_RUNS (subjectcode, session, output_dir)
%
% This function reads in the merged nifti of <subjectcode> <session> and 
% averages the timeseries of the two concatinated runs. The no. of dynamics
% of the averaged timeseries matches the run with the fewest dynamics. 
%
% Input: <subjectcode>, <session>, <output_dir>
% Output: <bairprf_AVERAGED_bold.nii>
%%

% ========= Read in nifti  ========= % 
merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_preproc.nii']);

run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold.nii']);
dyn_run1 = size(run1,4);
run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold.nii']);
dyn_run2 = size(run2,4);    
n_vol = min(dyn_run1, dyn_run2);

part1 = merged (:,:,:, 1:n_vol);
part2 = merged (:,:,:, (1+n_vol):(2*n_vol));
  
% ========= Take average of 2 pRF runs  ========= % 
avrg = (part1+part2)/2;

% ========= Write output nifti (task-bairprf_AVERAGED_preproc)  ========= % 
hdr = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold.nii']);
filename = [subjectcode, '_', session, '_task-bairprf_AVERAGED_preproc.nii'];
niftiwrite (avrg, fullfile(output_dir, filename), hdr)

%%

disp (['Done averaging pRF runs for ', subjectcode, ': ', session])
