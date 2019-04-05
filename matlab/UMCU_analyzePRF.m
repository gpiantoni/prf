%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%

%% Dataset: Visual01 3T

% Read in niftis (preprocessed) 
img = 'sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_bold-rwm.nii';
img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz';
run1 = niftiread (img);
img = '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz';
run2 = niftiread (img);

% Speficy vector of voxel indices to analyze
vxs = zeros(size(run1, 1), size(run1, 2), size(run1, 3));
vxs(30:60, 30:60, 26) = run1(30:60, 30:60, 26) >= 0;

% Load in apertures
apertures = load('/home/giovanni/tools/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');
images = uint8(apertures.bar_apertures); clear apertures;

temp = zeros(100,100,size(images,3));
for p=1:size(images,3)
    temp(:,:,p) = imresize(images(:,:,p),[100 100],'cubic');
end
images = temp;

% ensure that all values are between 0 and 1
images(images < 0) = 0;
images(images > 1) = 1;

% add baseline (first 11.9 sec)
TR = .850;
BASELINE = 11.9; % seconds
BASELINE_TR = 11.9 / TR;
stimulus_baseline = zeros(100, 100, BASELINE_TR);

images = cat(3, stimulus_baseline, images);
images = images(:, :, 1:224);

%%

% Compute pRF parameters with GLMDenoise
results  = analyzePRF({images, images},{run1, run2}, TR, struct('seedmode',[0 1],'display','on', 'wantglmdenoise', 1, 'vxs', find(vxs)));

% GLMDenoise returns noise regressors in <noisereg>

%%
info = niftiinfo(img);

info.ImageSize = info.ImageSize(1:3);
info.PixelDimensions = info.PixelDimensions(1:3);
info.Datatype = 'double';
niftiwrite(results.R2, '/home/margriet/Documents/MATLAB/data_kay/analyzePRF-master/analyzePRF_results/r2.nii', info);
niftiwrite(results.ang, '/home/margriet/Documents/MATLAB/data_kay/analyzePRF-master/analyzePRF_results/ang.nii', info);

load ('/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01.mat')
load ('/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/sub-visual01_ses-UMCU7TGE_task-bairtemporalpatterns_run-01.mat')