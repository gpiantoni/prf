% % =========  SCRIPT response buttonbox ========= % 
%
%%
clear; 

subjectcode = 'sub-visual12';

run1_3TMB = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01.mat']);
run2_3TMB = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02.mat']);
run1_7TGE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01.mat']);
run2_7TGE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02.mat']);
run1_7TSE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01.mat']);
run2_7TSE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02.mat']);

per_run1_3TMB = run1_3TMB.pc;
per_run2_3TMB = run2_3TMB.pc;
per_run1_7TGE = run1_7TGE.pc;
per_run2_7TGE = run2_7TGE.pc;
per_run1_7TSE = run1_7TSE.pc;
per_run2_7TSE = run2_7TSE.pc;

response = [per_run1_3TMB, per_run2_3TMB; per_run1_7TGE, per_run2_7TGE; per_run1_7TSE, per_run2_7TSE];

figure(1)
bar (response)
set (gca, 'XTickLabel', {'3TMB', '7TGE', '7TSE'})
ylabel ('Percentage')
title (['Response: ', subjectcode])
set (gcf, 'Position', [800, 800, 1000, 800])

mean = mean(response(:))

%% Save plots
% 
% outdir = '/Fridge/users/margriet/projects/prf/buttonbox_response';
% saveas (figure(1), fullfile(outdir, subjectcode), 'png')

%%

disp ('End')
