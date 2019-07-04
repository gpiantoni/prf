function convert_pixels2degrees (session, output_dir)
%
% function CONVERT_PIXELS2DEGREES (subject, session, output_dir)
%
% This function converts pixels to visual angle in degrees 
%
% Input: <subject>, <session>, <output_dir>
% Output: <ecc_deg> and <rfsize_deg> found in <output_dir>
%
%%
disp('Converting pixels to degrees visual angle')

%% Calculate visual angle

% ========= 3T scanner ========= % 
scanner3T = load ('/Fridge/users/margriet/stimuli/BAIR_pRF/scanner_params/Scanner_3T.mat');        % 3T scanner parameters
scanner3T.params.display.numPixels = [1200, 1920];          % pixels
scanner3T.params.display.dimensions = [32.40, 51.84];       % [32.40, 51.84] cm 
scanner3T.params.display.distance = 141.8;                  % 141.8 cm
scanner3T.params.display.pixelSize = scanner3T.params.display.dimensions / scanner3T.params.display.numPixels;  % 1 pixel = 0.027 X 0.027 cm
screen_size_3T = scanner3T.params.display.dimensions;
distance_3T = scanner3T.params.display.distance;

visual_angle_3T = (atan((screen_size_3T(1)/2) / distance_3T)) * 180/pi;      % in degrees (6.5175 deg)  5.59

% ========= 7T scanner ========= % 
scanner7T = load ('/Fridge/users/margriet/stimuli/BAIR_pRF/scanner_params/Scanner_7T.mat');        % 7T scanner parameters
screen_size_7T = scanner7T.params.display.dimensions;     % [14.22, 8] cm
screen_size_7T = [screen_size_7T(2), screen_size_7T(1)];  % [8, 14.22] cm
distance_7T = scanner7T.params.display.distance;          % 35.5 cm

visual_angle_7T = (atan((screen_size_7T(1)/2) / distance_7T)) * 180/pi;      % in degrees (6.428 deg)
   
%% Read in results (eccentricity and rf size)

ecc = niftiread (fullfile(output_dir, 'ecc.nii'));
rfsize = niftiread (fullfile(output_dir, 'rfsize.nii'));

%% Convert pixels to degrees

img_resolution = [100, 100];
pix_div2 = img_resolution(1)/2;

if session == 'ses-UMCU3TMB'
    ecc_deg = ecc * (visual_angle_3T/pix_div2);
    rfsize_deg = rfsize * (visual_angle_3T/pix_div2);
else
    ecc_deg = ecc * (visual_angle_7T/pix_div2);
    rfsize_deg = rfsize * (visual_angle_7T/pix_div2); 
end

%% Output

% ========= Read in header information ========= % 
hdr_ecc = niftiinfo (fullfile(output_dir, 'ecc.nii'));
hdr_rfsize = niftiinfo (fullfile(output_dir, 'rfsize.nii'));

% ========= Write output nifti file ========= % 
niftiwrite (ecc_deg, fullfile(output_dir, 'ecc_deg.nii'), hdr_ecc);
niftiwrite (rfsize_deg, fullfile(output_dir, 'rfsize_deg.nii'), hdr_rfsize);

%% End
