% % =========  SCRIPT pRF results by area  ========= % 
%
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

%% Specify parameters

subjectcode = 'sub-visual03'; 
subjectnumber = str2num(subjectcode (11:12));
session = 'ses-UMCU3TMB';

%%

% ========= READ IN PRF RESULTS ========= %
ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/ang_masked1.nii']);
ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/ecc_masked1.nii']);
rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/rfsize_masked1.nii']);
R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/R2.nii']);

% ========= READ IN BENSON ATLAS ========= %
varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_MERGED_bold/varea.nii']);
Benson_ROI_Names = {'V1', 'V2', 'V3', 'hV4', 'VO1', 'VO2', 'LO1', 'LO2', 'TO1', 'TO2', 'V3B', 'V3A'};

% ========= RESHAPE ========= %
dim = size (R2);       % [fx X fy X fz X dynamics]
fx = dim(1);
fy = dim(2);
fz = dim(3);
nrnodes = fx*fy*fz;          % To convert 4D to 2D = 276480 (3T)

ang_2d = reshape (ang, [nrnodes,1]);
ecc_2d = reshape (ecc, [nrnodes,1]);
rfsize_2d = reshape (rfsize, [nrnodes,1]);
varea_2d = reshape (varea, [nrnodes,1]);


% ========= RESULTS BY AREA ========= %
%%% POLAR ANGLE %%%
ang_V1 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 1
        ang_V1(i) = ang_2d(i);
    end
end

ang_V2 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 2
        ang_V2(i) = ang_2d(i);
    end
end

ang_V3 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 3
        ang_V3(i) = ang_2d(i);
    end
end

%%% ECCENTRICITY %%%
ecc_V1 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 1
        ecc_V1(i) = ecc_2d(i);
    end
end

ecc_V2 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 2
        ecc_V2(i) = ecc_2d(i);
    end
end

ecc_V3 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 3
        ecc_V3(i) = ecc_2d(i);
    end
end

%%% RECEPTIVE FIELD SIZE %%%
rfsize_V1 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 1
        rfsize_V1(i) = rfsize_2d(i);
    end
end

rfsize_V2 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 2
        rfsize_V2(i) = rfsize_2d(i);
    end
end

rfsize_V3 = zeros (nrnodes, 1);
for i = 1:nrnodes
    if varea_2d(i) == 3
        rfsize_V3(i) = rfsize_2d(i);
    end
end


% ========= PLOT RESULTS: SCATTERPLOTS ========= %

%%


