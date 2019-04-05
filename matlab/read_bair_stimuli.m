function images = read_bair_stimuli(n_volumes, TR)
%READ_BAIR_STIMULI read in the apertures with the correct TR and number of
%scans
%
% <n_volumes>


PATH_TO_APERTURES = '/home/giovanni/tools/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat';
BASELINE = 11.9; % seconds
RESOLUTION = 100;

apertures = load(PATH_TO_APERTURES);
images = uint8(apertures.bar_apertures);

temp = zeros(RESOLUTION, RESOLUTION, size(images,3));
for p=1:size(images,3)
    temp(:,:,p) = imresize(images(:,:,p), [RESOLUTION RESOLUTION], 'cubic');
end
images = temp;

% ensure that all values are between 0 and 1
images(images < 0) = 0;
images(images > 1) = 1;

% add baseline
BASELINE_TR = (BASELINE / TR);
stimulus_baseline = zeros(RESOLUTION, RESOLUTION, BASELINE_TR);

images = cat(3, stimulus_baseline, images);
images = images(:, :, 1:n_volumes);
