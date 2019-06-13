% % =========  SCRIPT pRF results by area  ========= % 
%
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

Analyze3TMB = 1;
Analyze7TGE = 0;
Analyze7TSE = 0;

Benson_ROI_Names = {'V1', 'V2', 'V3', 'hV4', 'VO1', 'VO2', 'LO1', 'LO2', 'TO1', 'TO2', 'V3B', 'V3A'};

%% Specify parameters

subjectcode = 'sub-visual03'; 
method = 'separate_bairprf';

V1_rfsize = (V1_rfsize - min(V1_rfsize)) / ( max(V1_rfsize) - min(V1_rfsize) );

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

    Results_3TMB = {};
    Results_3TMB.ang = ang_1d;
    Results_3TMB.ecc = ecc_1d;
    Results_3TMB.rfsize = rfsize_1d;
    Results_3TMB.R2 = R2_1d;

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/varea.nii']);
    varea_3TMB = reshape (varea, [nrnodes,1]);

    % ========= COMPUTE RESULTS BY VISUAL AREA ========= %
         %%% POLAR ANGLE %%%
    V1_ang = zeros (nrnodes, 1);            % V1
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_ang(i) = Results_3TMB.ang (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);            % V2
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_ang(i) = Results_3TMB.ang (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);            % V3
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_ang(i) = Results_3TMB.ang (i);
        end
    end

         %%% ECCENTRICITY %%%
    V1_ecc = zeros (nrnodes, 1);            % V1
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_ecc(i) = Results_3TMB.ecc(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);            % V2
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_ecc(i) = Results_3TMB.ecc(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);            % V3
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_ecc(i) = Results_3TMB.ecc(i);
        end
    end

         %%% RECEPTIVE FIELD SIZE %%%
    V1_rfsize = zeros (nrnodes, 1);         % V1
    for i = 1:nrnodes
        if varea_3TMB(i) == 1
            V1_rfsize(i) = Results_3TMB.rfsize(i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);         % V2
    for i = 1:nrnodes
        if varea_3TMB(i) == 2
            V2_rfsize(i) = Results_3TMB.rfsize(i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         % V3
    for i = 1:nrnodes
        if varea_3TMB(i) == 3
            V3_rfsize(i) = Results_3TMB.rfsize(i);
        end
    end
    
    % ========= REMOVE ZEROES ========= %
    V1_ang = V1_ang(V1_ang>0);
    V1_ecc = V1_ecc(V1_ecc>0);
    V1_rfsize = V1_rfsize(V1_rfsize>0);
    
     % ========= REMOVE INF ========= %
    V1_ang = V1_ang(V1_ang~=Inf);
    V1_ecc = V1_ecc(V1_ecc~=Inf);
    V1_rfsize = V1_rfsize(V1_rfsize~=Inf);
    
    
    
    V1_ang(V1_ang == Inf) = 999;
    V1_rfsize(V1_rfsize == Inf) = 999
    
    
    i = find (V1_rfsize>0);
    j = find (i(V1_rfsize(i) ~=Inf));
    k = find (j(~isnan(V1_rfsize(j))));
    
    
    a = find(isinf(V1_rfsize));
    b = find(isnan(V1_rfsize));
    
    index = cat (1, a, b);
    
    V1_rfsize(index) = [];
    V1_ecc (index) = [];
    V1_ang (index) = [];
   
    
    
    
    
%     k = find (V1_ang>0);
%     l = find (k(V1_ang(k) ~=Inf));
%     
%     m = find (V1_ecc>0);
%     n = find (m(V1_ecc(m) ~=Inf));
%     
    
    
    C = j ~= l
    
    
     
    275823
    
    A = [1 2 3 4 Inf 2 1];
    B = [1 2 3 4 7 2 1];
    
    C = A ~= B
    
    
       
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
    
    % ========= PLOT RESULTS ========= %
    P = polyfit(V1_3TMB.ecc, V1_3TMB.rfsize, 1);
    yfit = polyval (P, V1_3TMB.ecc);
    
    figure(1)
    plot (V1_3TMB.ecc, V1_3TMB.rfsize, '*r')
    hold on
    plot (V1_3TMB.ecc, yfit, '-k')
      xlabel ('Eccentricity')
    ylabel ('Receptive field size')
    title ('V1 (ses-UMCU3TMB)')
    axis([0 16 0 12])

    %%%%% V2 %%%%%
    figure(2)
    plot (V2_3TMB.ecc, V2_3TMB.rfsize, '*b')
    xlabel ('Eccentricity')
    ylabel ('Receptive field size')
    title ('V2 (ses-UMCU3TMB)')

    %%%%% V3 %%%%%
    figure(3)
    plot (V3_3TMB.ecc, V3_3TMB.rfsize, '*m')
    xlabel ('Eccentricity')
    ylabel ('Receptive field size')
    title ('V3 (ses-UMCU3TMB)')  
    
    
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

    Results_7TGE = {};
    Results_7TGE.ang = ang_1d;
    Results_7TGE.ecc = ecc_1d;
    Results_7TGE.rfsize = rfsize_1d;
    Results_7TGE.R2 = R2_1d;

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/varea.nii']);
    varea_7TGE = reshape (varea, [nrnodes,1]);
    
    % ========= RESULTS BY VISUAL AREA ========= %
         %%% POLAR ANGLE %%%
    V1_ang = zeros (nrnodes, 1);            % V1
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_ang(i) = Results_7TGE.ang (i);
        end
    end
    V2_ang = zeros (nrnodes, 1);            % V2
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_ang(i) = Results_7TGE.ang (i);
        end
    end
    V3_ang = zeros (nrnodes, 1);            % V3
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_ang(i) = Results_7TGE.ang (i);
        end
    end

         %%% ECCENTRICITY %%%
    V1_ecc = zeros (nrnodes, 1);            % V1
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_ecc(i) = Results_7TGE.ecc(i);
        end
    end
    V2_ecc = zeros (nrnodes, 1);            % V2
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_ecc(i) = Results_7TGE.ecc(i);
        end
    end
    V3_ecc = zeros (nrnodes, 1);            % V3
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_ecc(i) = Results_7TGE.ecc(i);
        end
    end

         %%% RECEPTIVE FIELD SIZE %%%
    V1_rfsize = zeros (nrnodes, 1);         % V1
    for i = 1:nrnodes
        if varea_7TGE(i) == 1
            V1_rfsize(i) = Results_7TGE.rfsize(i);
        end
    end
    V2_rfsize = zeros (nrnodes, 1);         % V2
    for i = 1:nrnodes
        if varea_7TGE(i) == 2
            V2_rfsize(i) = Results_7TGE.rfsize(i);
        end
    end
    V3_rfsize = zeros (nrnodes, 1);         % V3
    for i = 1:nrnodes
        if varea_7TGE(i) == 3
            V3_rfsize(i) = Results_7TGE.rfsize(i);
        end
    end
    
    % ========= PLOT RESULTS ========= %
    
    
    
    
    
    
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

    Results_7TSE = {};
    Results_7TSE.ang = ang_1d;
    Results_7TSE.ecc = ecc_1d;
    Results_7TSE.rfsize = rfsize_1d;
    Results_7TSE.R2 = R2_1d;

    % ========= READ IN BENSON ATLAS ========= %
    varea = niftiread (['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/varea.nii']);
    varea_7TSE = reshape (varea, [nrnodes,1]);
end

% % % %% Plotting data
% % % 
% % % % ========= RESULTS BY AREA ========= %
% % % %%% POLAR ANGLE %%%
% % % V1_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % V2_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % V3_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 3
% % %         V3_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % %%% ECCENTRICITY %%%
% % % V1_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % V2_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % V3_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 3
% % %         V3_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % %%% RECEPTIVE FIELD SIZE %%%
% % % V1_rfsize = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_rfsize(i) = rfsize_1d(i);
% % %     end
% % % end
% % % 
% % % V2_rfsize = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_rfsize(i) = rfsize_1d(i);
% % %     end
% % % end
% % % % ========= RESULTS BY AREA ========= %
% % % %%% POLAR ANGLE %%%
% % % V1_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % V2_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % V3_ang = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 3
% % %         V3_ang(i) = ang_1d(i);
% % %     end
% % % end
% % % 
% % % %%% ECCENTRICITY %%%
% % % V1_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % V2_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % V3_ecc = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 3
% % %         V3_ecc(i) = ecc_1d(i);
% % %     end
% % % end
% % % 
% % % %%% RECEPTIVE FIELD SIZE %%%
% % % V1_rfsize = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 1
% % %         V1_rfsize(i) = rfsize_1d(i);
% % %     end
% % % end
% % % 
% % % V2_rfsize = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 2
% % %         V2_rfsize(i) = rfsize_1d(i);
% % %     end
% % % end
% % % 
% % % V3_rfsize = zeros (nrnodes, 1);
% % % for i = 1:nrnodes
% % %     if varea_1d(i) == 3
% % %         V3_rfsize(i) = rfsize_1d(i);
% % %     end
% % % end

% ========= PLOT RESULTS: SCATTERPLOTS ========= %

%% Clear variables in workspace

% clear V1_ang
% clear V1_ecc
% clear V1_rfsize
% clear V2_ang
% clear V2_ecc
% clear V2_rfsize
% clear V3_ang
% clear V3_ecc
% clear V3_rfsize


