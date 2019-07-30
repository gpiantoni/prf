% % =========  SCRIPT results for R ========= % 
%
% Script to convert pRF results to csv files for analysis in R.
%
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

subjectcode = 'sub-visual11'; 
subjectname = subjectcode (5:end);
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
    V1_R2 = zeros (nrnodes, 1);         %%% R2 %%%
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_R2(i) = R2_1d (i);
        end
    end
    V2_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_R2(i) = R2_1d (i);
        end
    end
    V3_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_R2(i) = R2_1d(i);
        end
    end

    %% REMOVE NaN & Inf & zeroes & outliers FROM DATA
    
    % ========= REMOVE BAD DATA FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    indexZero = find(V1_ang == 0);
    indexOutlier = find (V1_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end
    
    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    indexZero = find(V1_ecc == 0);
    indexOutlier = find (V1_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
          
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    indexZero = find (V1_rfsize == 0);
    indexOutlier = find (V1_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
           
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
        
    %%% REMOVE SELECTED INDICES %%%
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];
    V1_R2 (rm_index) = [];
    
    % ========= REMOVE BAD DATA FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    indexZero = find(V2_ang == 0);
    indexOutlier = find (V2_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    indexZero = find(V2_ecc == 0);
    indexOutlier = find (V2_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end

    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    indexZero = find (V2_rfsize == 0);
    indexOutlier = find (V2_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);

    %%% REMOVE SELECTED INDICES %%%
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    V2_R2 (rm_index) = [];
    
    % ========= REMOVE BAD DATA FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    indexZero = find (V3_ang == 0);
    indexOutlier = find (V3_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    indexZero = find (V3_ecc == 0);
    indexOutlier = find (V3_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    indexZero = find (V3_rfsize == 0);
    indexOutlier = find (V3_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
    
    %%% REMOVE SELECTED INDICES %%%
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];     
    V3_R2 (rm_index) = [];
       
    %% GROUP RESULTS - 3TMB
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_3TMB = [V1_ecc, V1_rfsize, V1_ang, V1_R2];
    V2_3TMB = [V2_ecc, V2_rfsize, V2_ang, V2_R2];
    V3_3TMB = [V3_ecc, V3_rfsize, V3_ang, V3_R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_3TMB_V1'];
    csvwrite (filename, V1_3TMB)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_3TMB_V2'];
    csvwrite (filename, V2_3TMB)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_3TMB_V3'];
    csvwrite (filename, V3_3TMB)
            
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
    V1_R2 = zeros (nrnodes, 1);              %%% R2 %%%
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_R2(i) = R2_1d(i);
        end
    end
    V2_R2 = zeros (nrnodes, 1);        
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_R2(i) = R2_1d(i);
        end
    end
    V3_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_R2(i) = R2_1d(i);
        end
    end

    %% REMOVE NaN & Inf & zeroes & outliers FROM DATA
    
    % ========= REMOVE BAD DATA FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    indexZero = find (V1_ang == 0);
    indexOutlier = find (V1_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end    

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    indexZero = find (V1_ecc == 0);
    indexOutlier = find (V1_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    indexZero = find (V1_rfsize == 0);
    indexOutlier = find (V1_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
            
    %%% REMOVE SELECTED INDICES %%%
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];    
    V1_R2 (rm_index) = [];  
    
    % ========= REMOVE BAD DATA FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    indexZero = find (V2_ang == 0);
    indexOutlier = find (V2_ang > 360);    
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    indexZero = find (V2_ecc == 0);
    indexOutlier = find (V2_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    indexZero = find (V2_rfsize == 0);
    indexOutlier = find (V2_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end

    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
        
    %%% REMOVE SELECTED INDICES %%%
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    V2_R2 (rm_index) = [];  
  
    % ========= REMOVE BAD DATA FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    indexZero = find (V3_ang == 0);
    indexOutlier = find (V3_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    indexZero = find (V3_ecc == 0);
     indexOutlier = find (V3_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    indexZero = find (V3_rfsize == 0);
    indexOutlier = find (V3_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
    
    %%% REMOVE SELECTED INDICES %%%   
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];       
    V3_R2 (rm_index) = [];      
   
    %% GROUP RESULTS - 7TGE
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_7TGE = [V1_ecc, V1_rfsize, V1_ang, V1_R2];
    V2_7TGE = [V2_ecc, V2_rfsize, V2_ang, V2_R2];
    V3_7TGE = [V3_ecc, V3_rfsize, V3_ang, V3_R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TGE_V1'];
    csvwrite (filename, V1_7TGE)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TGE_V2'];
    csvwrite (filename, V2_7TGE)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TGE_V3'];
    csvwrite (filename, V3_7TGE)
    
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
    V1_R2 = zeros (nrnodes, 1);         %%% R2 %%%
    for i = 1:nrnodes
        if varea_7TSE(i) == 1
            V1_R2(i) = R2_1d(i);
        end
    end
    V2_R2 = zeros (nrnodes, 1);        
    for i = 1:nrnodes
        if varea_7TSE(i) == 2
            V2_R2(i) = R2_1d(i);
        end
    end
    V3_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea_7TSE(i) == 3
            V3_R2(i) = R2_1d(i);
        end
    end

    %% REMOVE NaN & Inf & zeroes & outliers FROM DATA
    
    % ========= REMOVE BAD DATA FOR V1 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V1_ang));
    indexNaN = find(isnan(V1_ang));
    indexZero = find (V1_ang == 0);
    indexOutlier = find (V1_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%%
    indexInf = find(isinf(V1_ecc));
    indexNaN = find(isnan(V1_ecc));
    indexZero = find (V1_ecc == 0);
    indexOutlier = find (V1_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%%
    indexInf = find(isinf(V1_rfsize));
    indexNaN = find(isnan(V1_rfsize));
    indexZero = find (V1_rfsize == 0);
    indexOutlier = find (V1_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
        
    %%% REMOVE SELECTED INDICES %%%   
    V1_rfsize(rm_index) = [];
    V1_ecc (rm_index) = [];
    V1_ang (rm_index) = [];
    V1_R2 (rm_index) = [];
    
    % ========= REMOVE BAD DATA FOR V2 ========= %
    
    %%% FIND BAD DATA IN ANG %%%
    indexInf = find(isinf(V2_ang));
    indexNaN = find(isnan(V2_ang));
    indexZero = find (V2_ang == 0);
    indexOutlier = find (V2_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V2_ecc));
    indexNaN = find(isnan(V2_ecc));
    indexZero = find (V2_ecc == 0);
    indexOutlier = find (V2_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V2_rfsize));
    indexNaN = find(isnan(V2_rfsize));
    indexZero = find (V2_rfsize == 0);
    indexOutlier = find (V2_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
        
    %%% REMOVE SELECTED INDICES %%%   
    V2_rfsize(rm_index) = [];
    V2_ecc (rm_index) = [];
    V2_ang (rm_index) = [];
    V2_R2 (rm_index) = [];
    
    % ========= REMOVE BAD DATA FOR V3 ========= %
    
    %%% FIND BAD DATA IN ANG %%
    indexInf = find(isinf(V3_ang));
    indexNaN = find(isnan(V3_ang));
    indexZero = find (V3_ang == 0);
    indexOutlier = find (V3_ang > 360);
    index_ang = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ang = indexInf;
    else 
        rm_Inf_ang = [];
    end

    %%% FIND BAD DATA IN ECC %%
    indexInf = find(isinf(V3_ecc));
    indexNaN = find(isnan(V3_ecc));
    indexZero = find (V3_ecc == 0);
    indexOutlier = find (V3_ecc > 25);
    index_ecc = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_ecc = indexInf;
    else 
        rm_Inf_ecc = [];
    end
    
    %%% FIND BAD DATA IN RFSIZE %%
    indexInf = find(isinf(V3_rfsize));
    indexNaN = find(isnan(V3_rfsize));
    indexZero = find (V3_rfsize == 0);
    indexOutlier = find (V3_rfsize > 20);
    index_rfsize = cat (1, indexInf, indexNaN, indexZero, indexOutlier);
    
    if ~isempty(indexInf)
        rm_Inf_rfsize = indexInf;
    else 
        rm_Inf_rfsize = [];
    end
    
    %%% SELECT LARGEST MATRIX TO REMOVE %%%
    indexdims = [length(index_ang), length(index_ecc), length(index_rfsize)];
    [maxdim, position] = max(indexdims);
    
    if position == 1 
        tmp_rm_index = index_ang;
    elseif position == 2
        tmp_rm_index = index_ecc;
    else 
        tmp_rm_index = index_rfsize;
    end
    
    %%% MAKE SURE INF VALUES ARE INCLUDED %%%
    tmp_rm_index_cat = cat (1, tmp_rm_index, rm_Inf_ang, rm_Inf_ecc, rm_Inf_rfsize);
    
    %%% REMOVE DUPLICATE VALUES %%%
    rm_index = unique (tmp_rm_index_cat);
        
    %%% REMOVE SELECTED INDICES %%%   
    V3_rfsize(rm_index) = [];
    V3_ecc (rm_index) = [];
    V3_ang (rm_index) = [];        
    V3_R2 (rm_index) = [];
     
       
    %% GROUP RESULTS - 7TSE
    
    % ========= RESULTS BY VISUAL AREA ========= %
    V1_7TSE = [V1_ecc, V1_rfsize, V1_ang, V1_R2];
    V2_7TSE = [V2_ecc, V2_rfsize, V2_ang, V2_R2];
    V3_7TSE = [V3_ecc, V3_rfsize, V3_ang, V3_R2];
    
    % SAVE RESULTS IN CSV
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TSE_V1'];
    csvwrite (filename, V1_7TSE)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TSE_V2'];
    csvwrite (filename, V2_7TSE)
    filename = ['/home/margriet/R/data/pRFresults/', subjectname, '_7TSE_V3'];
    csvwrite (filename, V3_7TSE)
    
end
