function split_merged_runs (subjectcode, session, output_dir)
%
% function SPLIT_MERGED_RUNS (subjectcode, session)
%
% This function reads in the merged nifti of <subjectcode> <session> and 
% splits up the concatinated runs to create 2 separate pRF runs. 
%
% Input: <subjectcode>, <session>, <output_dir>
% Output: <bairprf_MERGED_runxx_bold.nii>
%%

% ========= Read in nifti  ========= % 
merged = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_preproc.nii']);

run1 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode, '_', session, '_task-bairprf_run-01_bold.nii']);
dyn_run1 = size(run1,4);
run2 = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode, '_', session, '_task-bairprf_run-02_bold.nii']);
dyn_run2 = size(run2,4);    

n_vol = min(dyn_run1, dyn_run2);

merged_part1 = merged (:,:,:, 1:dyn_run1);
merged_part2 = merged (:,:,:, (1+dyn_run1):(dyn_run1+dyn_run2));

% ========= Write output nifti (task-bairprf_MERGED_runxx_preproc)  ========= % 
% Run1
hdr1 = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_preproc.nii']);
hdr1.ImageSize(4) = dyn_run1;
filename_run1 = [subjectcode, '_', session, '_task-bairprf_MERGED_run01_preproc.nii'];
niftiwrite (merged_part1, fullfile(output_dir, filename_run1), hdr1)

% Run2
hdr2 = niftiinfo (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session '_task-bairprf_MERGED_bold/', subjectcode, '_', session, '_task-bairprf_MERGED_preproc.nii']);
hdr2.ImageSize(4) = dyn_run2;
filename_run2 = [subjectcode, '_', session, '_task-bairprf_MERGED_run02_preproc.nii'];
niftiwrite (merged_part2, fullfile(output_dir, filename_run2), hdr2)



%%

% figure(1)
% plot (squeeze(merged_part1(71,69,1, :)), 'LineWidth', 2, 'Color', 'r')
% hold on 
% plot (squeeze(merged_part2(71,69,1, :)), 'LineWidth', 2, 'Color', 'm')
% plot (squeeze(merged(71,69,1, :)), 'LineWidth', 2, 'Color', 'c')
% legend ('Run1', 'Run2', 'Merged')
% % axis ([0 500 -20 20])
% hold off


disp (['Done splitting merged pRF runs for ', subjectcode, ': ', session])














