function [V1, V2, V3] = compute_pRFresults_byArea (varea, ang, ecc, rfsize, R2, nrnodes)
%
% function [V1, V2, V3] = COMPUTE_PRFRESULTS_BYAREA (varea, ang, ecc, rfsize, R2, nrnodes)
%
% This function groups the pRF results (<ang>, <ecc>, <rfsize>, <R2>) based 
% on  visual area (<varea>) as defined by the Benson atlas (Benson et al.,
% 2014). The function returns three structure arrays containing the pRF
% properties of the voxels in V1, V2 and V3, respectively.
%
% Input: <varea>, <ang>, <ecc>, <rfsize>, <R2>, <nrnodes>
% Output: <[V1, V2, V3]>
%%

    V1_ang = zeros (nrnodes, 1);           %%% POLAR ANGLE %%%
    for i = 1:nrnodes
        if varea(i) == 1
            V1_ang(i) = ang (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);            
    for i = 1:nrnodes
        if varea(i) == 2
            V2_ang(i) = ang (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea(i) == 3
            V3_ang(i) = ang (i);
        end
    end
    V1_ecc = zeros (nrnodes, 1);             %%% ECCENTRICITY %%%
    for i = 1:nrnodes
        if varea(i) == 1
            V1_ecc(i) = ecc(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea(i) == 2
            V2_ecc(i) = ecc(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);           
    for i = 1:nrnodes
        if varea(i) == 3
            V3_ecc(i) = ecc(i);
        end
    end
    V1_rfsize = zeros (nrnodes, 1);         %%% RECEPTIVE FIELD SIZE %%%
    for i = 1:nrnodes
        if varea(i) == 1
            V1_rfsize(i) = rfsize (i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea(i) == 2
            V2_rfsize(i) = rfsize (i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea(i) == 3
            V3_rfsize(i) = rfsize(i);
        end
    end
    V1_R2 = zeros (nrnodes, 1);         %%% R2 %%%
    for i = 1:nrnodes
        if varea(i) == 1
            V1_R2(i) = R2 (i);
        end
    end
    V2_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea(i) == 2
            V2_R2(i) = R2 (i);
        end
    end
    V3_R2 = zeros (nrnodes, 1);         
    for i = 1:nrnodes
        if varea(i) == 3
            V3_R2(i) = R2(i);
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
    
    V1 = {};
    V1.rfsize = V1_rfsize;
    V1.ecc = V1_ecc;
    V1.ang = V1_ang;
    V1.R2 = V1_R2;
    
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
        
    V2 = {};
    V2.rfsize = V2_rfsize;
    V2.ecc = V2_ecc;
    V2.ang = V2_ang;
    V2.R2 = V2_R2;
    
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
        
    V3 = {};
    V3.rfsize = V3_rfsize;
    V3.ecc = V3_ecc;
    V3.ang = V3_ang;
    V3.R2 = V3_R2;
    
    

