function convert_pixels2degrees (subject, session, output_dir)
%
% function CONVERT_PIXELS2DEGREES (subject, session, output_dir)
%
% This function converts pixels to visual angle in degrees 
%
% Input: <subject>, <session>, <output_dir>
% Output: <ecc_deg.nii> and <rfsize_deg> found in <output_dir>
%
%%

disp('Converting pixels to degrees visual angle')

%% Calculate visual angle

% ========= 3T scanner ========= % 
scanner3T = load ('/Fridge/R01_BAIR/visual_fmri/data/raw/visual01/3T_MB/log/sub-visual01_ses-umc3t02_task-prf_run-1.mat');          % 3T scanner parameters
screen_size_3T = scanner3T.stimulus.display.dimensions;     % [32.5, 52] cm
distance_3T = scanner3T.stimulus.display.distance;          % 112 cm

visual_angle_3T = (atan((screen_size_3T(1)/2) / distance_3T)) * 180/pi;      % in degrees (8.2555 deg)

% ========= 7T scanner ========= % 
scanner7T = load ('/Fridge/R01_BAIR/visual_fmri/data/raw/visual01/7T_GE/log/sub-visual01GE_ses-umc7t01_task-prf_run-1.mat');        % 7T scanner parameters
screen_size_7T = scanner7T.stimulus.display.dimensions;     % [14.22, 8] cm
screen_size_7T = [screen_size_7T(2), screen_size_7T(1)];    % [8, 14.22] cm
distance_7T = scanner7T.stimulus.display.distance;          % 35.5 cm

visual_angle_7T = (atan((screen_size_7T(1)/2) / distance_7T)) * 180/pi;      % in degrees (6.428 deg)
   
%% Read in results (eccentricity and rf size)

ecc = niftiread (['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subject, '/', session, '/merged_bairprf/ecc.nii']);
rfsize = niftiread (['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subject, '/', session, '/merged_bairprf/rfsize.nii']);


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
hdr_ecc = niftiinfo (['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subject, '/', session, '/merged_bairprf/ecc.nii']);
hdr_rfsize = niftiinfo (['/Fridge/users/margriet/projects/analysis/analyzeprf/results/umcu/', subject, '/', session, '/merged_bairprf/rfsize.nii']);

% ========= Write output nifti file ========= % 
niftiwrite (ecc_deg, fullfile(output_dir, 'ecc_deg.nii'), hdr_ecc);
niftiwrite (rfsize_deg, fullfile(output_dir, 'rfsize_deg.nii'), hdr_rfsize);

%% End
