function images = read_bair_stimuli(n_volumes, TR)
% 
% function images = READ_BAIR_STIMULI (n_volumes, TR)
%
% This function reads in the apertures of the stimuli presented, along with
% the the number of dynamics (<n_volumes>) and the repetition time (<TR>).
%
% Output: 
% <images> containing pRF stimulus information with added baseline.

PATH_TO_APERTURES = '/Fridge/users/margriet/stimuli/BAIR_pRF/bar_apertures.mat';
BASELINE = 11.9;     % seconds
RESOLUTION = 100;

apertures = load(PATH_TO_APERTURES);
images = uint8(apertures.bar_apertures);

temp = zeros(RESOLUTION, RESOLUTION, size(images,3));
for p=1:size(images,3)
    temp(:,:,p) = imresize(images(:,:,p), [RESOLUTION RESOLUTION], 'cubic');
end
images = temp;

% Ensure that all values are between 0 and 1
images(images < 0) = 0;
images(images > 1) = 1;

% Add baseline
BASELINE_TR = round(BASELINE / TR);
stimulus_baseline = zeros(RESOLUTION, RESOLUTION, BASELINE_TR);

images = cat(3, stimulus_baseline, images);
images = images(:, :, 1:n_volumes);
