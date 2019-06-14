% % =========  SCRIPT pRF results by area  ========= % 
%
% Script visualizes the rfsize/ecc relation in V1, V2 and V3.
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

subjectcode = 'sub-visual03'; 
method = 'separate_bairprf';

%% 3TMB DATA

if Analyze3TMB == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/ang_masked1.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/ecc_masked1.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/rfsize_masked1.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/R2.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/varea.nii']);
    varea_3TMB = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
         
    V1_ang = zeros (nrnodes, 1);           %%% POLAR ANGLE %%%
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_ang(i) = ang_1d (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_ang(i) = ang_1d (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_ang(i) = ang_1d (i);
        end
    end
    V1_ecc = zeros (nrnodes, 1);             %%% ECCENTRICITY %%%
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_ecc(i) = ecc_1d(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_ecc(i) = ecc_1d(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_ecc(i) = ecc_1d(i);
        end
    end
    V1_rfsize = zeros (nrnodes, 1);         %%% RECEPTIVE FIELD SIZE %%%
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_rfsize(i) = rfsize_1d (i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_rfsize(i) = rfsize_1d (i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_rfsize(i) = rfsize_1d(i);
        end
    end

    %% REMOVE NaN & Inf FROM DATA
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    %%% REMOVE SELECTED INDICES %%
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end

    %%% REMOVE SELECTED INDICES %%
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end

    %%% REMOVE SELECTED INDICES %%
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];             
       
    %% CREATE CELLS
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_3TMB = {};
    V1_3TMB.ang = V1_ang;
    V1_3TMB.ecc = V1_ecc;
    V1_3TMB.rfsize = V1_rfsize;
    
    V2_3TMB = {};
    V2_3TMB.ang = V2_ang;
    V2_3TMB.ecc = V2_ecc;
    V2_3TMB.rfsize = V2_rfsize;
    
    V3_3TMB = {};
    V3_3TMB.ang = V3_ang;
    V3_3TMB.ecc = V3_ecc;
    V3_3TMB.rfsize = V3_rfsize;
    
    %% VISUALIZE RESULTS - ses-UMCU3TMB
    
     % Linear regression
    P = polyfit(V1_3TMB.ecc, V1_3TMB.rfsize, 1);
    yfit_P = polyval (P, V1_3TMB.ecc);
    Q = polyfit(V2_3TMB.ecc, V2_3TMB.rfsize, 1);
    yfit_Q = polyval (Q, V2_3TMB.ecc);
    R = polyfit(V3_3TMB.ecc, V3_3TMB.rfsize, 1);
    yfit_R = polyval (R, V3_3TMB.ecc);
        
    % ========= SCATTERPLOT 3T MULTIBAND ========= %
    figure (1) 
        subplot (1,3,1)
    plot (V1_3TMB.ecc, V1_3TMB.rfsize, '*r')
    hold on
    plot (V1_3TMB.ecc, yfit_P, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V1 (ses-UMCU3TMB)')
    axis([0 16 0 12])
        subplot(1,3,2)
    plot (V2_3TMB.ecc, V2_3TMB.rfsize, '*b')
    hold on
    plot (V2_3TMB.ecc, yfit_Q, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU3TMB)')
    axis([0 16 0 12])
        subplot(1,3,3)
    plot (V3_3TMB.ecc, V3_3TMB.rfsize, '*m')
    hold on
    plot (V3_3TMB.ecc, yfit_R, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU3TMB)')  
    axis([0 16 0 12])
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
    axis([0 16 0 8])
    legend ('V1', 'V2', 'V3')
    set (gcf, 'Position', [800, 800, 600, 800])      
            
end

%% 7TGE DATA

if Analyze7TGE == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/ang_masked1.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/ecc_masked1.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/rfsize_masked1.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/R2.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/varea.nii']);
    varea_7TGE = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
         
    V1_ang = zeros (nrnodes, 1);            %%% POLAR ANGLE %%%
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_ang(i) = ang_1d (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_ang(i) = ang_1d (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_ang(i) = ang_1d (i);
        end
    end
    V1_ecc = zeros (nrnodes, 1);             %%% ECCENTRICITY %%%
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_ecc(i) = ecc_1d(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_ecc(i) = ecc_1d(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);          
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_ecc(i) = ecc_1d(i);
        end
    end
    V1_rfsize = zeros (nrnodes, 1);         %%% RECEPTIVE FIELD SIZE %%%
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_rfsize(i) = rfsize_1d(i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);        
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_rfsize(i) = rfsize_1d(i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_rfsize(i) = rfsize_1d(i);
        end
    end

    %% REMOVE NaN & Inf FROM DATA
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    %%% REMOVE SELECTED INDICES %%
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    %%% REMOVE SELECTED INDICES %%
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
    
    %%% REMOVE SELECTED INDICES %%    
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];        
     
       
    %% CREATE CELLS
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_7TGE = {};
    V1_7TGE.ang = V1_ang;
    V1_7TGE.ecc = V1_ecc;
    V1_7TGE.rfsize = V1_rfsize;
    
    V2_7TGE = {};
    V2_7TGE.ang = V2_ang;
    V2_7TGE.ecc = V2_ecc;
    V2_7TGE.rfsize = V2_rfsize;
    
    V3_7TGE = {};
    V3_7TGE.ang = V3_ang;
    V3_7TGE.ecc = V3_ecc;
    V3_7TGE.rfsize = V3_rfsize;
    
    %% VISUALIZE RESULTS - ses-UMCU7TGE

     % Linear regression
    S = polyfit(V1_7TGE.ecc, V1_7TGE.rfsize, 1);
    yfit_S = polyval (S, V1_7TGE.ecc);
    T = polyfit(V2_7TGE.ecc, V2_7TGE.rfsize, 1);
    yfit_T = polyval (T, V2_7TGE.ecc);
    U = polyfit(V3_7TGE.ecc, V3_7TGE.rfsize, 1);
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
    axis([0 16 0 12])
        subplot(1,3,2)
    plot (V2_7TGE.ecc, V2_7TGE.rfsize, '*b')
    hold on
    plot (V2_7TGE.ecc, yfit_T, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TGE)')
    axis([0 16 0 12])
        subplot(1,3,3)
    plot (V3_7TGE.ecc, V3_7TGE.rfsize, '*m')
    hold on
    plot (V3_7TGE.ecc, yfit_U, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU7TGE)')  
    axis([0 16 0 12])
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
    axis([0 16 0 8])
    legend ('V1', 'V2', 'V3')          
    set (gcf, 'Position', [800, 800, 600, 800])
    
end

%% 7TSE DATA

if Analyze7TSE == true
    % ========= READ IN PRF RESULTS ========= %
    ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/ang_masked1.nii']);
    ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/ecc_masked1.nii']);
    rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/rfsize_masked1.nii']);
    R2 = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/R2.nii']);

    % ========= RESHAPE ========= %
    nrnodes = size(R2,1) * size(R2,2) * size(R2,3);

    ang_1d = reshape (ang, [nrnodes,1]);
    ecc_1d = reshape (ecc, [nrnodes,1]);
    rfsize_1d = reshape (rfsize, [nrnodes,1]);
    R2_1d = reshape (R2, [nrnodes,1]);

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/varea.nii']);
    varea_7TSE = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
         
    V1_ang = zeros (nrnodes, 1);            %%% POLAR ANGLE %%%
    for i = 1:nrnodes
        if varea_7TSE(i) == 1
            V1_ang(i) = ang_1d (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_7TSE(i) == 2
            V2_ang(i) = ang_1d (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_7TSE(i) == 3
            V3_ang(i) = ang_1d (i);
        end
    end
    V1_ecc = zeros (nrnodes, 1);            %%% ECCENTRICITY %%%
    for i = 1:nrnodes
        if varea_7TSE(i) == 1
            V1_ecc(i) = ecc_1d(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea_7TSE(i) == 2
            V2_ecc(i) = ecc_1d(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea_7TSE(i) == 3
            V3_ecc(i) = ecc_1d(i);
        end
    end
    V1_rfsize = zeros (nrnodes, 1);         %%% RECEPTIVE FIELD SIZE %%%
    for i = 1:nrnodes
        if varea_7TSE(i) == 1
            V1_rfsize(i) = rfsize_1d(i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);        
    for i = 1:nrnodes
        if varea_7TSE(i) == 2
            V2_rfsize(i) = rfsize_1d(i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_7TSE(i) == 3
            V3_rfsize(i) = rfsize_1d(i);
        end
    end

    %% REMOVE NaN & Inf FROM DATA
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    
    % ========= REMOVE BAD DATA (NaN & Inf) FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    index_ang = cat (1, indexInf, indexNaN);

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    index_ecc = cat (1, indexInf, indexNaN);
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    index_rfsize = cat (1, indexInf, indexNaN);
    
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim position] = max(indexdims);
    
    if position == 1 
        rm_index = index_ang;
    elseif position == 2
        rm_index = index_ecc;
    else 
        rm_index = index_rfsize;
    end
        
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];        
     
       
    %% CREATE CELLS
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_7TSE = {};
    V1_7TSE.ang = V1_ang;
    V1_7TSE.ecc = V1_ecc;
    V1_7TSE.rfsize = V1_rfsize;
    
    V2_7TSE = {};
    V2_7TSE.ang = V2_ang;
    V2_7TSE.ecc = V2_ecc;
    V2_7TSE.rfsize = V2_rfsize;
    
    V3_7TSE = {};
    V3_7TSE.ang = V3_ang;
    V3_7TSE.ecc = V3_ecc;
    V3_7TSE.rfsize = V3_rfsize;
    
    %% VISUALIZE RESULTS - ses-UMCU7TSE
    
    % Linear regression
    V = polyfit(V1_7TSE.ecc, V1_7TSE.rfsize, 1);
    yfit_V = polyval (V, V1_7TSE.ecc);
    W = polyfit(V2_7TSE.ecc, V2_7TSE.rfsize, 1);
    yfit_W = polyval (W, V2_7TSE.ecc);
    X = polyfit(V3_7TSE.ecc, V3_7TSE.rfsize, 1);
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
    axis([0 16 0 12])
        subplot(1,3,2)
    plot (V2_7TSE.ecc, V2_7TSE.rfsize, '*b')
    hold on
    plot (V2_7TSE.ecc, yfit_W, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V2 (ses-UMCU7TSE)')
    axis([0 16 0 12])
        subplot(1,3,3)
    plot (V3_7TSE.ecc, V3_7TSE.rfsize, '*m')
    hold on
    plot (V3_7TSE.ecc, yfit_X, '-k')
    hold off
    xlabel ('Eccentricity (^{o})')
    ylabel ('Receptive field size (^{o})')
    title ('V3 (ses-UMCU7TSE)')  
    axis([0 16 0 12])
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
    axis([0 16 0 8])
    legend ('V1', 'V2', 'V3')
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
axis([0 16 0 8])
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
axis([0 16 0 8])
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
axis([0 16 0 8])
legend ('ses-UMCU3TMB', 'ses-UMCU7TGE', 'ses-UMCU7TSE', 'Location', 'Northwest') 
set (gcf, 'Position', [800, 800, 600, 800])


%% Save plots



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
legend ('V1', 'V2', 'V3')

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
ylabel ('Standard deviation receptive field size (^{o})')
title (' receptive field size per visual area')  
legend ('V1', 'V2', 'V3')


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
clear index_ang
clear index_ecc
clear index_rfsize
clear indexdims
clear indexInf
clear indexNaN
clear maxdim
clear position
clear varea
clear i
clear mean*
clear std*

%% END