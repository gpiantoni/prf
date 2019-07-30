% % =========  SCRIPT pRF results by area  ========= % 
%
% Script visualizes the rfsize/ecc relation in V1, V2 and V3 based on the 
% pRF results thresholded with 1% R2.
%
% Input: <subjectcode>, <method>, <Analyze3TMB>, <Analyze7TGE>, <Analyze7TSE>
% Output: rfsize/ecc relation plots per area & rfsize histograms per area

%%
clear;
close all

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

Analyze3TMB = 1;
Analyze7TGE = 1;
Analyze7TSE = 1;

Benson_ROI_Names = {'V1', 'V2', 'V3', 'hV4', 'VO1', 'VO2', 'LO1', 'LO2', 'TO1', 'TO2', 'V3B', 'V3A'};

%% Specify parameters

subjectcode = 'sub-visual01'; 
method = 'final';

disp ('Running')

%% Color definition

color_area = [
    0.929   0.694   0.125;  % V1
    0.850   0.325   0.098;  % V2
    0.635   0.078   0.184]; % V3

color_seq = [
    0.3467    0.5360    0.6907      % 3TMB
    0.9153    0.2816    0.2878      % 7TGE
    0.4416    0.7490    0.4322];    % 7TSE

%% 3TMB DATA

if Analyze3TMB == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_3TMB/ang_FOVmasked.nii.gz']);
    ecc = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_3TMB/ecc_FOVmasked.nii.gz']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_3TMB/rfsize_FOVmasked.nii.gz']);
    R2 = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_3TMB/R2_FOVmasked.nii.gz']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);
    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/Benson/', subjectcode, '/BENSON-wide_3TMB_func.nii.gz']);
    varea_3TMB = reshape (varea, [nrnodes,1]);
    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
    [V1_3TMB, V2_3TMB, V3_3TMB] = compute_pRFresults_byArea (varea_3TMB, ang_1d, ecc_1d, rfsize_1d, R2_1d, nrnodes); 
    
    %% VISUALIZE RESULTS - ses-UMCU3TMB
    
    % Linear regression - Weighted R2
    P = polyfitweighted(V1_3TMB.ecc, V1_3TMB.rfsize, 1, V1_3TMB.R2);
    yfit_P = polyval (P, V1_3TMB.ecc);
    Q = polyfitweighted(V2_3TMB.ecc, V2_3TMB.rfsize, 1, V2_3TMB.R2);
    yfit_Q = polyval (Q, V2_3TMB.ecc);
    R = polyfitweighted(V3_3TMB.ecc, V3_3TMB.rfsize, 1, V3_3TMB.R2);
    yfit_R = polyval (R, V3_3TMB.ecc);
              
    % ========= SCATTERPLOT 3T MULTIBAND ========= %
   
    figure (1) 
        subplot (1,3,1)
    plot (V1_3TMB.ecc, V1_3TMB.rfsize, '*', 'Color', color_area(1,:))
    hold on
    plot (V1_3TMB.ecc, yfit_P, '-k')
%    errorbar(V1_3TMB.BINS_CENTER, V1_3TMB.rfsize_mean, V1_3TMB.rfsize_sem, 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.3, 'Color', 'k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU3TMB)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_3TMB.ecc, V2_3TMB.rfsize, '*', 'Color', color_area(2,:))
    hold on
    plot (V2_3TMB.ecc, yfit_Q, '-k')
%    errorbar(V2_3TMB.BINS_CENTER, V2_3TMB.rfsize_mean, V2_3TMB.rfsize_sem, 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.3, 'Color', 'k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU3TMB)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_3TMB.ecc, V3_3TMB.rfsize, '*', 'Color', color_area(3,:))
    hold on
    plot (V3_3TMB.ecc, yfit_R, '-k')
%    errorbar(V3_3TMB.BINS_CENTER, V3_3TMB.rfsize_mean, V3_3TMB.rfsize_sem, 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.3, 'Color', 'k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU3TMB)')  
    axis([0 12 0 6])
    set (gcf, 'Position', [800, 800, 1900, 400])
        
    % ========= SUMMARY PLOT - 3T MULTIBAND ========= %
    figure (2)
    plot (V1_3TMB.ecc, yfit_P, 'Color', color_area(1,:))
    hold on
    plot (V2_3TMB.ecc, yfit_Q, 'Color', color_area(2,:))
    plot (V3_3TMB.ecc, yfit_R, 'Color', color_area(3,:))
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ([subjectcode, ' (ses-UMCU3TMB)'])  
    axis([0 12 0 6])
    legend ('V1', 'V2', 'V3', 'Location', 'Northwest') 
    set (gcf, 'Position', [800, 800, 600, 800])   
    
          
end

%% 7TGE DATA

if Analyze7TGE == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TGE/ang_FOVmasked.nii.gz']);
    ecc = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TGE/ecc_FOVmasked.nii.gz']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TGE/rfsize_FOVmasked.nii.gz']);
    R2 = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TGE/R2_FOVmasked.nii.gz']);
    
    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);
    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/Benson/', subjectcode, '/BENSON-wide_7TGE_func.nii.gz']);
    varea_7TGE = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
    [V1_7TGE, V2_7TGE, V3_7TGE] = compute_pRFresults_byArea (varea_7TGE, ang_1d, ecc_1d, rfsize_1d, R2_1d, nrnodes);
    
        
    %% BIN RESULTS - ses-UMCU7TGE
   
%     [V1_7TGE.BINS_CENTER, V1_7TGE.rfsize_mean, V1_7TGE.rfsize_sem] = compute_bins(V1_7TGE.ecc, V1_7TGE.rfsize);
%     [V2_7TGE.BINS_CENTER, V2_7TGE.rfsize_mean, V2_7TGE.rfsize_sem] = compute_bins(V2_7TGE.ecc, V2_7TGE.rfsize);   
%     [V3_7TGE.BINS_CENTER, V3_7TGE.rfsize_mean, V3_7TGE.rfsize_sem] = compute_bins(V3_7TGE.ecc, V3_7TGE.rfsize);  
        
    %% VISUALIZE RESULTS - ses-UMCU7TGE

%      % Linear regression
%     S = polyfit(V1_7TGE.ecc, V1_7TGE.rfsize, 1);
%     yfit_S = polyval (S, V1_7TGE.ecc);
%     T = polyfit(V2_7TGE.ecc, V2_7TGE.rfsize, 1);
%     yfit_T = polyval (T, V2_7TGE.ecc);
%     U = polyfit(V3_7TGE.ecc, V3_7TGE.rfsize, 1);
%     yfit_U = polyval (U, V3_7TGE.ecc);
    
    % Linear regression - Weighted R2
    S = polyfitweighted(V1_7TGE.ecc, V1_7TGE.rfsize, 1, V1_7TGE.R2);
    yfit_S = polyval (S, V1_7TGE.ecc);
    T = polyfitweighted(V2_7TGE.ecc, V2_7TGE.rfsize, 1, V2_7TGE.R2);
    yfit_T = polyval (T, V2_7TGE.ecc);
    U = polyfitweighted(V3_7TGE.ecc, V3_7TGE.rfsize, 1, V3_7TGE.R2);
    yfit_U = polyval (U, V3_7TGE.ecc);
 
    % ========= SCATTERPLOT 7T GRADIENT ECHO ========= %
    figure (3) 
        subplot (1,3,1)
    plot (V1_7TGE.ecc, V1_7TGE.rfsize, '*', 'Color', color_area(1,:))
    hold on
    plot (V1_7TGE.ecc, yfit_S, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU7TGE)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_7TGE.ecc, V2_7TGE.rfsize, '*', 'Color', color_area(2,:))
    hold on
    plot (V2_7TGE.ecc, yfit_T, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TGE)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_7TGE.ecc, V3_7TGE.rfsize, '*', 'Color', color_area(3,:))
    hold on
    plot (V3_7TGE.ecc, yfit_U, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU7TGE)')  
    axis([0 12 0 6])
    set (gcf, 'Position', [800, 800, 1900, 400])
        
    % ========= SUMMARY PLOT - 7T GRADIENT ECHO ========= %
    figure (4)
    plot (V1_7TGE.ecc, yfit_S, 'Color', color_area(1,:))
    hold on
    plot (V2_7TGE.ecc, yfit_T, 'Color', color_area(2,:))
    plot (V3_7TGE.ecc, yfit_U, 'Color', color_area(3,:))
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ([subjectcode, ' (ses-UMCU7TGE)'])  
    axis([0 12 0 6])
    legend ('V1', 'V2', 'V3', 'Location', 'Northwest')           
    set (gcf, 'Position', [800, 800, 600, 800])
    
end

%% 7TSE DATA

if Analyze7TSE == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TSE/ang_FOVmasked.nii.gz']);
    ecc = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TSE/ecc_FOVmasked.nii.gz']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TSE/rfsize_FOVmasked.nii.gz']);
    R2 = niftiread (['/Fridge/users/margriet/projects/FOV/', subjectcode, '/pRF_7TSE/R2_FOVmasked.nii.gz']);
    
    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/Benson/', subjectcode, '/BENSON-wide_7TSE_func.nii.gz']);
    varea_7TSE = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
    [V1_7TSE, V2_7TSE, V3_7TSE] = compute_pRFresults_byArea (varea_7TSE, ang_1d, ecc_1d, rfsize_1d, R2_1d, nrnodes);
    
      %% BIN RESULTS - ses-UMCU7TSE
   
%     [V1_7TSE.BINS_CENTER, V1_7TSE.rfsize_mean, V1_7TSE.rfsize_sem] = compute_bins(V1_7TSE.ecc, V1_7TSE.rfsize);
%     [V2_7TSE.BINS_CENTER, V2_7TSE.rfsize_mean, V2_7TSE.rfsize_sem] = compute_bins(V2_7TSE.ecc, V2_7TSE.rfsize);   
%     [V3_7TSE.BINS_CENTER, V3_7TSE.rfsize_mean, V3_7TSE.rfsize_sem] = compute_bins(V3_7TSE.ecc, V3_7TSE.rfsize);    
        
    %% VISUALIZE RESULTS - ses-UMCU7TSE
    
%     % Linear regression
%     V = polyfit(V1_7TSE.ecc, V1_7TSE.rfsize, 1);
%     yfit_V = polyval (V, V1_7TSE.ecc);
%     W = polyfit(V2_7TSE.ecc, V2_7TSE.rfsize, 1);
%     yfit_W = polyval (W, V2_7TSE.ecc);
%     X = polyfit(V3_7TSE.ecc, V3_7TSE.rfsize, 1);
%     yfit_X = polyval (X, V3_7TSE.ecc);
    
    % Linear regression - Weighted R2
    V = polyfitweighted(V1_7TSE.ecc, V1_7TSE.rfsize, 1, V1_7TSE.R2);
    yfit_V = polyval (V, V1_7TSE.ecc);
    W = polyfitweighted(V2_7TSE.ecc, V2_7TSE.rfsize, 1, V2_7TSE.R2);
    yfit_W = polyval (W, V2_7TSE.ecc);
    X = polyfitweighted(V3_7TSE.ecc, V3_7TSE.rfsize, 1, V3_7TSE.R2);
    yfit_X = polyval (X, V3_7TSE.ecc);
    
    % ========= SCATTERPLOT - 7T SPIN ECHO ========= %
    figure (5) 
        subplot (1,3,1)
    plot (V1_7TSE.ecc, V1_7TSE.rfsize, '*', 'Color', color_area(1,:))
    hold on
    plot (V1_7TSE.ecc, yfit_V, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU7TSE)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_7TSE.ecc, V2_7TSE.rfsize, '*', 'Color', color_area(2,:))
    hold on
    plot (V2_7TSE.ecc, yfit_W, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TSE)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_7TSE.ecc, V3_7TSE.rfsize, '*', 'Color', color_area(3,:))
    hold on
    plot (V3_7TSE.ecc, yfit_X, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU7TSE)')  
    axis([0 12 0 6])
    set (gcf, 'Position', [800, 800, 1900, 400])
        
    % ========= SUMMARY PLOT - 7T SPIN ECHO ========= %
    figure (6)
    plot (V1_7TSE.ecc, yfit_V, 'Color', color_area(1,:))
    hold on
    plot (V2_7TSE.ecc, yfit_W, 'Color', color_area(2,:))
    plot (V3_7TSE.ecc, yfit_X, 'Color', color_area(3,:))
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ([subjectcode, ' (ses-UMCU7TSE)'])  
    axis([0 12 0 6])
    legend ('V1', 'V2', 'V3', 'Location', 'Northwest')
    set (gcf, 'Position', [800, 800, 600, 800])
     
end

%% Plot results for all sessions

% ========= SUMMARY PLOT - V1 ========= %
figure (7)
plot (V1_3TMB.ecc, yfit_P,  'Color', color_seq(1,:))
hold on
plot (V1_7TGE.ecc, yfit_S,  'Color', color_seq(2,:))
plot (V1_7TSE.ecc, yfit_V,  'Color', color_seq(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V1)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])

% ========= SUMMARY PLOT - V2 ========= %
figure (8)
plot (V2_3TMB.ecc, yfit_Q, 'Color', color_seq(1,:))
hold on
plot (V2_7TGE.ecc, yfit_T, 'Color', color_seq(2,:))
plot (V2_7TSE.ecc, yfit_W, 'Color', color_seq(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V2)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])

% ========= SUMMARY PLOT - V3 ========= %
figure (9)
plot (V3_3TMB.ecc, yfit_R, 'Color', color_seq(1,:))
hold on
plot (V3_7TGE.ecc, yfit_U, 'Color', color_seq(2,:))
plot (V3_7TSE.ecc, yfit_X, 'Color', color_seq(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V3)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])

%% PREP RESULTS FOR R 

% ========= 3TMB data for R ========= %
    % ========= RESULTS BY VISUAL AREA ========= %
    Rdata_3TMB_V1 = [V1_3TMB.ecc, V1_3TMB.rfsize, V1_3TMB.ang, V1_3TMB.R2];
    Rdata_3TMB_V2 = [V2_3TMB.ecc, V2_3TMB.rfsize, V2_3TMB.ang, V2_3TMB.R2];
    Rdata_3TMB_V3 = [V3_3TMB.ecc, V3_3TMB.rfsize, V3_3TMB.ang, V3_3TMB.R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V1_3TMB'];
    csvwrite (filename, Rdata_3TMB_V1)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V2_3TMB'];
    csvwrite (filename, Rdata_3TMB_V2)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V3_3TMB'];
    csvwrite (filename, Rdata_3TMB_V3)
    
% ========= 7TGE data for R ========= %
    % ========= RESULTS BY VISUAL AREA ========= %
    Rdata_7TGE_V1 = [V1_7TGE.ecc, V1_7TGE.rfsize, V1_7TGE.ang, V1_7TGE.R2];
    Rdata_7TGE_V2 = [V2_7TGE.ecc, V2_7TGE.rfsize, V2_7TGE.ang, V2_7TGE.R2];
    Rdata_7TGE_V3 = [V3_7TGE.ecc, V3_7TGE.rfsize, V3_7TGE.ang, V3_7TGE.R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V1_7TGE'];
    csvwrite (filename, Rdata_7TGE_V1)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V2_7TGE'];
    csvwrite (filename, Rdata_7TGE_V2)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V3_7TGE'];
    csvwrite (filename, Rdata_7TGE_V3)
    
    
% ========= 7TSE data for R ========= %
    % ========= RESULTS BY VISUAL AREA ========= %
    Rdata_7TSE_V1 = [V1_7TSE.ecc, V1_7TSE.rfsize, V1_7TSE.ang, V1_7TSE.R2];
    Rdata_7TSE_V2 = [V2_7TSE.ecc, V2_7TSE.rfsize, V2_7TSE.ang, V2_7TSE.R2];
    Rdata_7TSE_V3 = [V3_7TSE.ecc, V3_7TSE.rfsize, V3_7TSE.ang, V3_7TSE.R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V1_7TSE'];
    csvwrite (filename, Rdata_7TSE_V1)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V2_7TSE'];
    csvwrite (filename, Rdata_7TSE_V2)
    filename = ['/home/margriet/R/data/pRF_adjFOV/', subjectcode, '/V3_7TSE'];
    csvwrite (filename, Rdata_7TSE_V3)



%% Clear variables in workspace

clear V1_ang
clear V1_ecc
clear V1_rfsize
clear V2_ang
clear V2_ecc
clear V2_rfsize
clear V3_ang
clear V3_ecc
clear V3_rfsize
clear ang
clear ecc
clear rfsize
clear ang_1d
clear ecc_1d
clear rfsize_1d
clear index*
clear maxdim
clear position
clear varea
clear i
clear mean*
clear std*
clear tmp*
clear rm*

%% END

