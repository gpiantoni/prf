% % =========  SCRIPT response buttonbox ========= % 
%
%%
clear; 

subjectcode = 'sub-visual02';

%%%% 3T DATA %%%%
HRF_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairhrfpattern_run-01.mat']);
PRF_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01.mat']);
PRF_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02.mat']);
SPAT_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialpatterns_run-01.mat']);
SPAT_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialpatterns_run-02.mat']);
TMP_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-01.mat']);
TMP_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-02.mat']);
TMP_3T_run3 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-03.mat']);
TMP_3T_run4 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-04.mat']);

pcHRF_3T_run1 = HRF_3T_run1.pc;
pcPRF_3T_run1 = PRF_3T_run1.pc;
pcPRF_3T_run2 = PRF_3T_run2.pc;
pcSPAT_3T_run1 = SPAT_3T_run1.pc;
pcSPAT_3T_run2 = SPAT_3T_run2.pc;
pcTMP_3T_run1 = TMP_3T_run1.pc;
pcTMP_3T_run2 = TMP_3T_run2.pc;
pcTMP_3T_run3 = TMP_3T_run3.pc;
pcTMP_3T_run4 = TMP_3T_run4.pc;

performance3T = [pcHRF_3T_run1, pcPRF_3T_run1, pcPRF_3T_run2, pcSPAT_3T_run1, pcSPAT_3T_run2, pcTMP_3T_run1, pcTMP_3T_run2, pcTMP_3T_run3, pcTMP_3T_run4];

%%%% 7TGE DATA %%%%
HRF_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairhrfpattern_run-01.mat']);
PRF_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01.mat']);
PRF_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02.mat']);
SPAT_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialpatterns_run-01.mat']);
SPAT_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialpatterns_run-02.mat']);
TMP_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-01.mat']);
TMP_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-02.mat']);
TMP_7TGE_run3 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-03.mat']);
TMP_7TGE_run4 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-04.mat']);

pcHRF_7TGE_run1 = HRF_7TGE_run1.pc;
pcPRF_7TGE_run1 = PRF_7TGE_run1.pc;
pcPRF_7TGE_run2 = PRF_7TGE_run2.pc;
pcSPAT_7TGE_run1 = SPAT_7TGE_run1.pc;
pcSPAT_7TGE_run2 = SPAT_7TGE_run2.pc;
pcTMP_7TGE_run1 = TMP_7TGE_run1.pc;
pcTMP_7TGE_run2 = TMP_7TGE_run2.pc;
pcTMP_7TGE_run3 = TMP_7TGE_run3.pc;
pcTMP_7TGE_run4 = TMP_7TGE_run4.pc;

performance7TGE = [pcHRF_7TGE_run1, pcPRF_7TGE_run1, pcPRF_7TGE_run2, pcSPAT_7TGE_run1, pcSPAT_7TGE_run2, pcTMP_7TGE_run1, pcTMP_7TGE_run2, pcTMP_7TGE_run3, pcTMP_7TGE_run4];

%%%% 7TSE DATA %%%%
HRF_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairhrfpattern_run-01.mat']);
PRF_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01.mat']);
PRF_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02.mat']);
SPAT_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialpatterns_run-01.mat']);
SPAT_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialpatterns_run-02.mat']);
TMP_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-01.mat']);
TMP_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-02.mat']);
TMP_7TSE_run3 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-03.mat']);
TMP_7TSE_run4 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-04.mat']);

pcHRF_7TSE_run1 = HRF_7TSE_run1.pc;
pcPRF_7TSE_run1 = PRF_7TSE_run1.pc;
pcPRF_7TSE_run2 = PRF_7TSE_run2.pc;
pcSPAT_7TSE_run1 = SPAT_7TSE_run1.pc;
pcSPAT_7TSE_run2 = SPAT_7TSE_run2.pc;
pcTMP_7TSE_run1 = TMP_7TSE_run1.pc;
pcTMP_7TSE_run2 = TMP_7TSE_run2.pc;
pcTMP_7TSE_run3 = TMP_7TSE_run3.pc;
pcTMP_7TSE_run4 = TMP_7TSE_run4.pc;

performance7TSE = [pcHRF_7TSE_run1, pcPRF_7TSE_run1, pcPRF_7TSE_run2, pcSPAT_7TSE_run1, pcSPAT_7TSE_run2, pcTMP_7TSE_run1, pcTMP_7TSE_run2, pcTMP_7TSE_run3, pcTMP_7TSE_run4];


%% SUBJECT VISUAL01
% % 
% % subjectcode = 'sub-visual01';
% % 
% % %%%% 3T DATA %%%%
% % HRF_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairhrfpattern_run-01.mat']);
% % PRF_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01.mat']);
% % PRF_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02.mat']);
% % SPAT_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialpatterns_run-01.mat']);
% % SPAT_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialpatterns_run-02.mat']);
% % TMP_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-01.mat']);
% % TMP_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairtemporalpatterns_run-02.mat']);
% % SPATOB_3T_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialobjects_run-01.mat']);
% % SPATOB_3T_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairspatialobjects_run-02.mat']);
% % 
% % pcHRF_3T_run1 = HRF_3T_run1.pc;
% % pcPRF_3T_run1 = PRF_3T_run1.pc;
% % pcPRF_3T_run2 = PRF_3T_run2.pc;
% % pcSPAT_3T_run1 = SPAT_3T_run1.pc;
% % pcSPAT_3T_run2 = SPAT_3T_run2.pc;
% % pcTMP_3T_run1 = TMP_3T_run1.pc;
% % pcTMP_3T_run2 = TMP_3T_run2.pc;
% % pcSPATOB_3T_run1 = SPATOB_3T_run1.pc;
% % pcSPATOB_3T_run2 = SPATOB_3T_run2.pc;
% % 
% % performance3T = [pcHRF_3T_run1, pcPRF_3T_run1, pcPRF_3T_run2, pcSPAT_3T_run1, pcSPAT_3T_run2, pcTMP_3T_run1, pcTMP_3T_run2, pcSPATOB_3T_run1, pcSPATOB_3T_run2];
% % 
% % %%%% 7TGE DATA %%%%
% % HRF_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairhrfpattern_run-01.mat']);
% % PRF_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01.mat']);
% % PRF_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02.mat']);
% % SPAT_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialpatterns_run-01.mat']);
% % SPAT_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialpatterns_run-02.mat']);
% % TMP_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-01.mat']);
% % TMP_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairtemporalpatterns_run-02.mat']);
% % SPATOB_7TGE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialobjects_run-01.mat']);
% % SPATOB_7TGE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairspatialobjects_run-02.mat']);
% % 
% % pcHRF_7TGE_run1 = HRF_7TGE_run1.pc;
% % pcPRF_7TGE_run1 = PRF_7TGE_run1.pc;
% % pcPRF_7TGE_run2 = PRF_7TGE_run2.pc;
% % pcSPAT_7TGE_run1 = SPAT_7TGE_run1.pc;
% % pcSPAT_7TGE_run2 = SPAT_7TGE_run2.pc;
% % pcTMP_7TGE_run1 = TMP_7TGE_run1.pc;
% % pcTMP_7TGE_run2 = TMP_7TGE_run2.pc;
% % pcSPATOB_7TGE_run1 = SPATOB_7TGE_run1.pc;
% % pcSPATOB_7TGE_run2 = SPATOB_7TGE_run2.pc;
% % 
% % performance7TGE = [pcHRF_7TGE_run1, pcPRF_7TGE_run1, pcPRF_7TGE_run2, pcSPAT_7TGE_run1, pcSPAT_7TGE_run2, pcTMP_7TGE_run1, pcTMP_7TGE_run2, pcSPATOB_7TGE_run1, pcSPATOB_7TGE_run2];
% % 
% % %%%% 7TSE DATA %%%%
% % HRF_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairhrfpattern_run-01.mat']);
% % PRF_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01.mat']);
% % PRF_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02.mat']);
% % SPAT_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialpatterns_run-01.mat']);
% % SPAT_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialpatterns_run-02.mat']);
% % TMP_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-01.mat']);
% % TMP_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairtemporalpatterns_run-02.mat']);
% % SPATOB_7TSE_run1 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialobjects_run-01.mat']);
% % SPATOB_7TSE_run2 = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairspatialobjects_run-02.mat']);
% % 
% % pcHRF_7TSE_run1 = HRF_7TSE_run1.pc;
% % pcPRF_7TSE_run1 = PRF_7TSE_run1.pc;
% % pcPRF_7TSE_run2 = PRF_7TSE_run2.pc;
% % pcSPAT_7TSE_run1 = SPAT_7TSE_run1.pc;
% % pcSPAT_7TSE_run2 = SPAT_7TSE_run2.pc;
% % pcTMP_7TSE_run1 = TMP_7TSE_run1.pc;
% % pcTMP_7TSE_run2 = TMP_7TSE_run2.pc;
% % pcSPATOB_7TSE_run1 = SPATOB_7TSE_run1.pc;
% % pcSPATOB_7TSE_run2 = SPATOB_7TSE_run2.pc;
% % 
% % performance7TSE = [pcHRF_7TSE_run1, pcPRF_7TSE_run1, pcPRF_7TSE_run2, pcSPAT_7TSE_run1, pcSPAT_7TSE_run2, pcTMP_7TSE_run1, pcTMP_7TSE_run2, pcSPATOB_7TSE_run1, pcSPATOB_7TSE_run2];

%% 
% % run1_3TMB = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-01.mat']);
% % run2_3TMB = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU3TMB_task-bairprf_run-02.mat']);
% % run1_7TGE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-01.mat']);
% % run2_7TGE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TGE_task-bairprf_run-02.mat']);
% % run1_7TSE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-01.mat']);
% % run2_7TSE = load (['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_ses-UMCU7TSE_task-bairprf_run-02.mat']);
% % 
% % per_run1_3TMB = run1_3TMB.pc;
% % per_run2_3TMB = run2_3TMB.pc;
% % per_run1_7TGE = run1_7TGE.pc;
% % per_run2_7TGE = run2_7TGE.pc;
% % per_run1_7TSE = run1_7TSE.pc;
% % per_run2_7TSE = run2_7TSE.pc;
% % 
% % response = [per_run1_3TMB, per_run2_3TMB; per_run1_7TGE, per_run2_7TGE; per_run1_7TSE, per_run2_7TSE];
% % 
% % figure(1)
% % bar (response)
% % set (gca, 'XTickLabel', {'3TMB', '7TGE', '7TSE'})
% % ylabel ('Percentage')
% % title (['Response: ', subjectcode])
% % set (gcf, 'Position', [800, 800, 1000, 800])
% % 
% % mean = mean(response(:))

%% Save plots
% 
% outdir = '/Fridge/users/margriet/projects/prf/buttonbox_response';
% saveas (figure(1), fullfile(outdir, subjectcode), 'png')

%%

disp ('End')
