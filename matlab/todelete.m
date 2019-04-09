% testscript

% subjectcode = 'sub-visual01'
% 
% 
% nifti =  {'/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz',
%     '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz'}
% 
% nifti2 =  {['/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz'],
%     ['/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz']}

subjectcode = 'sub-visual02';
number = str2num(subjectcode (11:12));

% Threshold NYU preprocessed data
thresholdNYU_3TMB = [1200];
thresholdNYU_7TGE = [];
thresholdNYU_7TSE = [];

% Threshold UMCU preprocessed data
thresholdUMCU_3TMB = [500, 750];
thresholdUMCU_7TGE = [150, 100];
thresholdUMCU_7TSE = [];

threshold = thresholdUMCU_3TMB(number);