function output_images = read_bair_stimuli(subjectcode, session, n_volumes, TR);
% 
% function images = READ_BAIR_STIMULI (n_volumes, TR)
%
% This function reads in the apertures of the stimuli presented, along with
% the the number of dynamics (<n_volumes>) and the repetition time (<TR>).
% <n_volumes> is a list of integers, with the number of volumes for each
% session
%
% Output: 
% <output_img> containing pRF stimulus information with added baseline.
%%

% % % % apertures = load(PATH_TO_APERTURES);
% % % % images = uint8(apertures.bar_apertures);

% % % % ========= Ensure that all images are binary ========= % 
% % % 
% % % images_bin = zeros(size(images));
% % % images_bin(images > .5) = 1;
% % % images = images_bin;


% n_volumes = 224;
% TR = 0.85;
% session = 'ses-UMCU7TGE';
% subjectcode = 'sub-visual02';

% ========= Read in stimuli per run ========= % 
PATH_TO_STIMULI_run1 = ['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_', session, '_task-bairprf_run-01.mat'];      % First run
img = load (PATH_TO_STIMULI_run1);
images_run1 = img.stimulus.images;
PATH_TO_STIMULI_run2 = ['/Fridge/R01_BAIR/visual_fmri/data/bids/stimuli/', subjectcode, '_', session, '_task-bairprf_run-02.mat'];      % Second run
img = load (PATH_TO_STIMULI_run2);
images_run2 = img.stimulus.images;

% ========= Turn stimuli into binary mask: include all non-grey pixels  ========= % 
images_run1 = images_run1 ~= 128;
images_run2 = images_run2 ~= 128;

% ========= Convert to 100X100 resolution ========= % 
RESOLUTION = 100;

% Run1
temp = zeros(RESOLUTION, RESOLUTION, size(images_run1,3));
for p=1:size(images_run1,3)
    temp(:,:,p) = imresize(images_run1(:,:,p), [RESOLUTION RESOLUTION], 'cubic');
end
images_run1 = temp;

% Run2
temp = zeros(RESOLUTION, RESOLUTION, size(images_run2,3));
for p=1:size(images_run2,3)
    temp(:,:,p) = imresize(images_run2(:,:,p), [RESOLUTION RESOLUTION], 'cubic');
end
images_run2 = temp;

% ========= Add baseline ========= % 
BASELINE = 11.9;                         % In seconds
BASELINE_TR = round(BASELINE / TR);      % In dynamics
stimulus_baseline = zeros(RESOLUTION, RESOLUTION, BASELINE_TR);

images_run1 = cat(3, stimulus_baseline, images_run1, stimulus_baseline);
images_run2 = cat(3, stimulus_baseline, images_run2, stimulus_baseline);


% ========= OUTPUT with n_volumes ========= % 
output_img_run1 = {};
for i = 1:length(n_volumes)
    output_img_run1{i} = images_run1(:, :, 1:n_volumes(i));
end
output_img_run1 = cat(3, output_img_run1{:});

output_img_run2 = {};
for i = 1:length(n_volumes)
    output_img_run2{i} = images_run2(:, :, 1:n_volumes(i));
end
output_img_run2 = cat(3, output_img_run2{:});

output_images = {output_img_run1, output_img_run2};

%% End