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

subjectcode = 'sub-visual12'; 
method = 'final';

disp ('Running')

%% 3TMB DATA

if Analyze3TMB == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/ang_masked5clust.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/ecc_masked5clust.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/rfsize_masked5clust.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/R2_masked5clust.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/anat2anat/', subjectcode, '/BENSON_varea_3TMB_func.nii.gz']);
    varea_3TMB = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
    [V1_3TMB, V2_3TMB, V3_3TMB] = compute_pRFresults_byArea (varea_3TMB, ang_1d, ecc_1d, rfsize_1d, R2_1d, nrnodes); 

    
    %% BIN RESULTS - ses-UMCU3TMB
   
%     [BINS_CENTER, rfsize_mean, rfsize_sem] = compute_bins(V1_3TMB.ecc, V1_3TMB.rfsize);
%     V1_3TMB.BINS_CENTER= BINS_CENTER;
%     V1_3TMB.rfsize_mean = rfsize_mean;
%     V1_3TMB.rfsize_sem = rfsize_sem;
%     
%     [V2_3TMB.BINS_CENTER, V2_3TMB.rfsize_mean, V2_3TMB.rfsize_sem] = compute_bins(V2_3TMB.ecc, V2_3TMB.rfsize);   
%     [V3_3TMB.BINS_CENTER, V3_3TMB.rfsize_mean, V3_3TMB.rfsize_sem] = compute_bins(V3_3TMB.ecc, V3_3TMB.rfsize);    
    
    %% VISUALIZE RESULTS - ses-UMCU3TMB
    
%      % Linear regression
%     P = polyfit(V1_3TMB.ecc, V1_3TMB.rfsize, 1);
%     yfit_P = polyval (P, V1_3TMB.ecc);
%     Q = polyfit(V2_3TMB.ecc, V2_3TMB.rfsize, 1);
%     yfit_Q = polyval (Q, V2_3TMB.ecc);
%     R = polyfit(V3_3TMB.ecc, V3_3TMB.rfsize, 1);
%     yfit_R = polyval (R, V3_3TMB.ecc);
    
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
    plot (V1_3TMB.ecc, V1_3TMB.rfsize, '*r')
    hold on
    plot (V1_3TMB.ecc, yfit_P, '-k')
%    errorbar(V1_3TMB.BINS_CENTER, V1_3TMB.rfsize_mean, V1_3TMB.rfsize_sem, 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.3, 'Color', 'k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU3TMB)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_3TMB.ecc, V2_3TMB.rfsize, '*b')
    hold on
    plot (V2_3TMB.ecc, yfit_Q, '-k')
%    errorbar(V2_3TMB.BINS_CENTER, V2_3TMB.rfsize_mean, V2_3TMB.rfsize_sem, 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.3, 'Color', 'k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU3TMB)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_3TMB.ecc, V3_3TMB.rfsize, '*m')
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
    plot (V1_3TMB.ecc, yfit_P, '-r')
    hold on
    plot (V2_3TMB.ecc, yfit_Q, '-b')
    plot (V3_3TMB.ecc, yfit_R, '-m')
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
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/ang_masked5clust.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/ecc_masked5clust.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/rfsize_masked5clust.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/R2_masked5clust.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/anat2anat/', subjectcode, '/BENSON_varea_7TGE_func.nii.gz']);
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
    plot (V1_7TGE.ecc, V1_7TGE.rfsize, '*r')
    hold on
    plot (V1_7TGE.ecc, yfit_S, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU7TGE)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_7TGE.ecc, V2_7TGE.rfsize, '*b')
    hold on
    plot (V2_7TGE.ecc, yfit_T, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TGE)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_7TGE.ecc, V3_7TGE.rfsize, '*m')
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
    plot (V1_7TGE.ecc, yfit_S, '-r')
    hold on
    plot (V2_7TGE.ecc, yfit_T, '-b')
    plot (V3_7TGE.ecc, yfit_U, '-m')
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
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/ang_masked5clust.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/ecc_masked5clust.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/rfsize_masked5clust.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/R2_masked5clust.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/projects/coregistration/anat2anat/', subjectcode, '/BENSON_varea_7TSE_func.nii.gz']);
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
    plot (V1_7TSE.ecc, V1_7TSE.rfsize, '*r')
    hold on
    plot (V1_7TSE.ecc, yfit_V, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU7TSE)')
    axis([0 12 0 6])
        subplot(1,3,2)
    plot (V2_7TSE.ecc, V2_7TSE.rfsize, '*b')
    hold on
    plot (V2_7TSE.ecc, yfit_W, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TSE)')
    axis([0 12 0 6])
        subplot(1,3,3)
    plot (V3_7TSE.ecc, V3_7TSE.rfsize, '*m')
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
    plot (V1_7TSE.ecc, yfit_V, '-r')
    hold on
    plot (V2_7TSE.ecc, yfit_W, '-b')
    plot (V3_7TSE.ecc, yfit_X, '-m')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ([subjectcode, ' (ses-UMCU7TSE)'])  
    axis([0 12 0 6])
    legend ('V1', 'V2', 'V3', 'Location', 'Northwest')
    set (gcf, 'Position', [800, 800, 600, 800])
     
end

%% Plot results for all sessions

color = [
    0.3467    0.5360    0.6907      % 3TMB
    0.9153    0.2816    0.2878      % 7TGE
    0.4416    0.7490    0.4322];    % 7TSE

% ========= SUMMARY PLOT - V1 ========= %
figure (7)
plot (V1_3TMB.ecc, yfit_P,  'Color', color(1,:))
hold on
plot (V1_7TGE.ecc, yfit_S,  'Color', color(2,:))
plot (V1_7TSE.ecc, yfit_V,  'Color', color(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V1)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])

% ========= SUMMARY PLOT - V2 ========= %
figure (8)
plot (V2_3TMB.ecc, yfit_Q, 'Color', color(1,:))
hold on
plot (V2_7TGE.ecc, yfit_T, 'Color', color(2,:))
plot (V2_7TSE.ecc, yfit_W, 'Color', color(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V2)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])

% ========= SUMMARY PLOT - V3 ========= %
figure (9)
plot (V3_3TMB.ecc, yfit_R, 'Color', color(1,:))
hold on
plot (V3_7TGE.ecc, yfit_U, 'Color', color(2,:))
plot (V3_7TSE.ecc, yfit_X, 'Color', color(3,:))
hold off
xlabel ('Eccentricity (^{o})')
ylabel ('Receptive field size (^{o})')
title ([subjectcode, ' (V3)'])  
axis([0 12 0 6])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])


%% Compute mean rfsize per sequence per area

% ========= MEAN RFSIZE - 3TMB ========= %
mean_3TMB_V1 = mean (V1_3TMB.rfsize);
mean_3TMB_V2 = mean (V2_3TMB.rfsize);
mean_3TMB_V3 = mean (V3_3TMB.rfsize);
MEAN_3TMB_rfsize = [mean_3TMB_V1, mean_3TMB_V2, mean_3TMB_V3];

% ========= MEAN RFSIZE - 7TGE ========= %
mean_7TGE_V1 = mean (V1_7TGE.rfsize);
mean_7TGE_V2 = mean (V2_7TGE.rfsize);
mean_7TGE_V3 = mean (V3_7TGE.rfsize);
MEAN_7TGE_rfsize = [mean_7TGE_V1, mean_7TGE_V2, mean_7TGE_V3];

% ========= MEAN RFSIZE - 7TSE ========= %
mean_7TSE_V1 = mean (V1_7TSE.rfsize);
mean_7TSE_V2 = mean (V2_7TSE.rfsize);
mean_7TSE_V3 = mean (V3_7TSE.rfsize);
MEAN_7TSE_rfsize = [mean_7TSE_V1, mean_7TSE_V2, mean_7TSE_V3];

MEAN_rfsize = [MEAN_3TMB_rfsize; MEAN_7TGE_rfsize; MEAN_7TSE_rfsize];

%%% Plot barcharts %%%
figure (10)
bar (MEAN_rfsize)
set (gca, 'XTickLabel', {'3TMB', '7TGE', '7TSE'})
ylabel ('Receptive field size (^{o})')
title ('Mean receptive field size per visual area')  

%% Compute standard deviation rfsize per sequence per area

% ========= STD RFSIZE - 3TMB ========= %
std_3TMB_V1 = std (V1_3TMB.rfsize);
std_3TMB_V2 = std (V2_3TMB.rfsize);
std_3TMB_V3 = std (V3_3TMB.rfsize);
STD_3TMB_rfsize = [std_3TMB_V1, std_3TMB_V2, std_3TMB_V3];

% ========= STD RFSIZE - 7TGE ========= %
std_7TGE_V1 = std (V1_7TGE.rfsize);
std_7TGE_V2 = std (V2_7TGE.rfsize);
std_7TGE_V3 = std (V3_7TGE.rfsize);
STD_7TGE_rfsize = [std_7TGE_V1, std_7TGE_V2, std_7TGE_V3];

% ========= STD RFSIZE - 7TSE ========= %
std_7TSE_V1 = std (V1_7TSE.rfsize);
std_7TSE_V2 = std (V2_7TSE.rfsize);
std_7TSE_V3 = std (V3_7TSE.rfsize);
STD_7TSE_rfsize = [std_7TSE_V1, std_7TSE_V2, std_7TSE_V3];

STD_rfsize = [STD_3TMB_rfsize; STD_7TGE_rfsize; STD_7TSE_rfsize];

%%% Plot barcharts %%%
figure (11)
bar (STD_rfsize)
set (gca, 'XTickLabel', {'3TMB', '7TGE', '7TSE'})
ylabel ('Receptive field size (^{o})')
title ('Standard deviation receptive field size per visual area')  

%% Plot histograms

% ========= HISTOGRAM - 3TMB - EXCL zeroes ========= %
figure (12)
subplot (1,3,1)
histogram (V1_3TMB.rfsize(V1_3TMB.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V1 (ses-UMCU3TMB)')
subplot (1,3,2)
histogram (V2_3TMB.rfsize(V2_3TMB.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V2 (ses-UMCU3TMB)')
subplot(1,3,3)
histogram (V3_3TMB.rfsize(V3_3TMB.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V3 (ses-UMCU3TMB)')

% ========= HISTOGRAM - 7TGE - EXCL zeroes ========= %
figure (13)
subplot (1,3,1)
histogram (V1_7TGE.rfsize(V1_7TGE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V1 (ses-UMCU7TGE)')
subplot (1,3,2)
histogram (V2_7TGE.rfsize(V2_7TGE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V2 (ses-UMCU7TGE)')
subplot(1,3,3)
histogram (V3_7TGE.rfsize(V3_7TGE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V3 (ses-UMCU7TGE)')

% ========= HISTOGRAM - 7TSE - EXCL zeroes ========= %
figure (14)
subplot (1,3,1)
histogram (V1_7TSE.rfsize(V1_7TSE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V1 (ses-UMCU7TSE)')
subplot (1,3,2)
histogram (V2_7TSE.rfsize(V2_7TSE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V2 (ses-UMCU7TSE)')
subplot(1,3,3)
histogram (V3_7TSE.rfsize(V3_7TSE.rfsize>0))
axis ([0 15 0 350])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title ('V3 (ses-UMCU7TSE)')

% ========= HISTOGRAM - SUMMARY ========= %
figure (15)
subplot (1,3,1)
histogram (V1_3TMB.rfsize(V1_3TMB.rfsize>0.5))
hold on
histogram (V1_7TGE.rfsize(V1_7TGE.rfsize>0.5))
histogram (V1_7TSE.rfsize(V1_7TSE.rfsize>0.5))
hold off
axis ([0 15 0 200])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title (['V1 (', subjectcode, ')'])
legend ('3TMB', '7TGE', '7TSE')
subplot (1,3,2)
histogram (V2_3TMB.rfsize(V2_3TMB.rfsize>0.5))
hold on
histogram (V2_7TGE.rfsize(V2_7TGE.rfsize>0.5))
histogram (V2_7TSE.rfsize(V2_7TSE.rfsize>0.5))
hold off
axis ([0 15 0 200])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title (['V2 (', subjectcode, ')'])
legend ('3TMB', '7TGE', '7TSE')
subplot (1,3,3)
histogram (V3_3TMB.rfsize(V3_3TMB.rfsize>0.5))
hold on
histogram (V3_7TGE.rfsize(V3_7TGE.rfsize>0.5))
histogram (V3_7TSE.rfsize(V3_7TSE.rfsize>0.5))
hold off
axis ([0 15 0 200])
xlabel ('Receptive field size (^{o})')
ylabel ('Frequency')
title (['V3 (', subjectcode, ')'])
legend ('3TMB', '7TGE', '7TSE')

%% Plot histogram eccentricity

% ========= HISTOGRAM - V1 eccentricity ========= %
figure (16)
a = histfit (V1_3TMB.ecc,50);
hold on
b = histfit (V1_7TGE.ecc,50);
c = histfit (V1_7TSE.ecc,50);
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8)
set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8)
set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8)
set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V1 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
set (gcf, 'Position', [800, 800, 600, 1200])

% ========= HISTOGRAM - V2 eccentricity ========= %
figure (17)
a = histfit (V2_3TMB.ecc, 50, 'normal');
hold on
b = histfit (V2_7TGE.ecc, 50, 'normal');
c = histfit (V2_7TSE.ecc, 50, 'normal');
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8)
set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8)
set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8)
set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V2 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
set (gcf, 'Position', [800, 800, 600, 1200])

% ========= HISTOGRAM - V3 eccentricity ========= %
figure (18)
a = histfit (V3_3TMB.ecc, 50, 'normal');
hold on
b = histfit (V3_7TGE.ecc, 50, 'normal');
c = histfit (V3_7TSE.ecc, 50, 'normal');
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8); set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8); set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8); set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V3 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
set (gcf, 'Position', [800, 800, 600, 1200])

% ========= HISTOGRAM - Eccentricity ========= %
figure(19)
subplot (1,3,1)
a = histfit (V1_3TMB.ecc,50);
hold on
b = histfit (V1_7TGE.ecc,50);
c = histfit (V1_7TSE.ecc,50);
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8); set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8); set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8); set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V1 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
subplot(1,3,2)
a = histfit (V2_3TMB.ecc, 50, 'normal');
hold on
b = histfit (V2_7TGE.ecc, 50, 'normal');
c = histfit (V2_7TSE.ecc, 50, 'normal');
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8); set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8); set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8); set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V2 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
subplot (1,3,3)
a = histfit (V3_3TMB.ecc, 50, 'normal');
hold on
b = histfit (V3_7TGE.ecc, 50, 'normal');
c = histfit (V3_7TSE.ecc, 50, 'normal');
hold off
axis ([0 15 0 100])
set(a(1), 'facecolor', color(1,:), 'FaceAlpha', 0.8); set(a(2), 'color', 'k')
set(b(1), 'facecolor', color(2,:), 'FaceAlpha', 0.8); set(b(2), 'color', 'k')
set(c(1), 'facecolor', color(3,:), 'FaceAlpha', 0.8); set(c(2), 'color', 'k')
xlabel ('Eccentricity (^{o})')
ylabel ('Frequency')
title ('V3 eccentricity distribution')
legend ([a(1), b(1), c(1)], '3TMB', '7TGE', '7TSE')
set (gcf, 'Position', [800, 800, 1800, 800])


%% Save plots

%% Save plots (figure 4 & 5)

% Define path
figpath = '/Fridge/users/margriet/projects/figures/byarea/';

% Define filename
filename1 = [subjectcode, '_3TMB_byarea'];
filename2 = [subjectcode, '_3TMB_byarea_summary'];
filename3 = [subjectcode, '_7TGE_byarea'];
filename4 = [subjectcode, '_7TGE_byarea_summary'];
filename5 = [subjectcode, '_7TSE_byarea'];
filename6 = [subjectcode, '_7TSE_byarea_summary'];
filename10 = [subjectcode, '_mean_byarea'];
filename11 = [subjectcode, '_std_byarea'];
filename19 = [subjectcode, '_eccdistribution'];

% Save plots
saveas (figure(1), fullfile(figpath, filename1), 'png');
saveas (figure(2), fullfile(figpath, filename2), 'png');
saveas (figure(3), fullfile(figpath, filename3), 'png');
saveas (figure(4), fullfile(figpath, filename4), 'png');
saveas (figure(5), fullfile(figpath, filename5), 'png');
saveas (figure(6), fullfile(figpath, filename6), 'png');
saveas (figure(10), fullfile(figpath, filename10), 'png');
saveas (figure(11), fullfile(figpath, filename11), 'png');
saveas (figure(19), fullfile(figpath, filename19), 'png');


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

