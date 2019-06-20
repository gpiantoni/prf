% % =========  SCRIPT GLMdenoise UMCU data  ========= % 
%
%
%%
clear;

addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf/scripts'))
addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf'))
addpath(genpath('/home/margriet/tools/prf/matlab'))  

stimdur = 0.5;
tr = 0.85;

%%

subjectcode = 'sub-visual01'; 
subjectnumber = str2num(subjectcode (11:12));
method = 'separate_bairprf';

Analyze3TMB = 1;
Analyze7TGE = 1;
Analyze7TSE = 1;

%%

if Analyze3TMB == true
    session = 'ses-UMCU3TMB';
    disp('%%%%%%%%%%%% Starting GLMdenoise for UMCU 3T (MB) %%%%%%%%%%%%')
    
    nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_run01_bold-rwm.nii'],
             ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU3TMB/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU3TMB_task-bairprf_MERGED_run02_bold-rwm.nii']};

    data = {};
    for i = 1:length(nifti)
        data{i} = niftiread (nifti{i});
    end
    
    tsv = {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-01_events.tsv'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-02_events.tsv']};
   
    % Run1 
    events1 = tdfread(tsv{1});
    nrscans1 = size(data{1}, 4);
    onset1 = round(events1.onset / tr);
    baseline1 = find(events1.trial_type == 255);
    
    desmat1 = zeros(nrscans1, 1);   
    desmat1(onset1) = 1;
    desmat1(baseline1) = 0;
    
    desmat1 = desmat1(1:nrscans1, 1);
    
    % Run2
    events2 = tdfread(tsv{2});
    nrscans2 = size(data{2}, 4);
    onset2 = round(events2.onset / tr);
    baseline2 = find(events2.trial_type == 255);
    
    desmat2 = zeros(nrscans2, 1);   
    desmat2(onset2) = 1;
    desmat2(baseline2) = 0;
    
    desmat2 = desmat2(1:nrscans2, 1);
    
    % Combine design matrices
    design = {};
    design = {desmat1, desmat2};
    
    figuredir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method, '/GLM_figures'];
    
    [results_GLMdenoise, data_GLMdenoise] = GLMdenoisedata(design,data,stimdur,tr, [], [] , [], figuredir);
    
    disp('%%%%%%%%%%%% Done: UMCU 3T (MB) %%%%%%%%%%%%')
    
    % Save whole results structure
    output_dir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method];
    cd (output_dir)
    save('data_preGLMdenoise.mat', 'data')    
    save('data_postGLMdenoise.mat', 'data_GLMdenoise')                             % Save results struct in output_dir
    save('results_GLMdenoise.mat', '-struct', 'results_GLMdenoise')         % Save results struct in output_dir
    cd ('/home/margriet/tools/prf/matlab')
           
end

if Analyze7TGE == true
    session = 'ses-UMCU7TGE';
    disp('%%%%%%%%%%%% Starting GLMdenoise for UMCU 7T (GE) %%%%%%%%%%%%')
        
    nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_run01_bold-masked-mc-warp.nii'],
             ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TGE/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TGE_task-bairprf_MERGED_run02_bold-masked-mc-warp.nii']};

    data = {};
    for i = 1:length(nifti)
        data{i} = niftiread (nifti{i});
    end
    
    tsv = {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-01_events.tsv'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-02_events.tsv']};
   
    % Run1 
    events1 = tdfread(tsv{1});
    nrscans1 = size(data{1}, 4);
    onset1 = round(events1.onset / tr);
    baseline1 = find(events1.trial_type == 255);
    
    desmat1 = zeros(nrscans1, 1);   
    desmat1(onset1) = 1;
    desmat1(baseline1) = 0;
    
    desmat1 = desmat1(1:nrscans1, 1);

    % Run2
    events2 = tdfread(tsv{2});
    nrscans2 = size(data{2}, 4);
    onset2 = round(events2.onset / tr);
    baseline2 = find(events2.trial_type == 255);
    
    desmat2 = zeros(nrscans2, 1);   
    desmat2(onset2) = 1;
    desmat2(baseline2) = 0;
    
    desmat2 = desmat2(1:nrscans2, 1);
       
    % Combine design matrices
    design = {};
    design = {desmat1, desmat2};
    
    figuredir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method, '/GLM_figures'];
    
    [results_GLMdenoise, data_GLMdenoise] = GLMdenoisedata(design,data,stimdur,tr, [], [] , [], figuredir);
    
    disp('%%%%%%%%%%%% Done: UMCU 7T (GE) %%%%%%%%%%%%')
    
    % Save whole results structure
    output_dir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method];
    cd (output_dir)
    save('data_preGLMdenoise.mat', 'data')    
    save('data_postGLMdenoise.mat', 'data_GLMdenoise')                             % Save results struct in output_dir
    save('results_GLMdenoise.mat', '-struct', 'results_GLMdenoise')            % Save results struct in output_dir
    cd ('/home/margriet/tools/prf/matlab')
         
end

if Analyze7TSE == true
    session = 'ses-UMCU7TSE';
    disp('%%%%%%%%%%%% Starting GLMdenoise for UMCU 7T (SE) %%%%%%%%%%%%')
        
    nifti =  {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_run01_bold-masked-mc-warp.nii'],
             ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/ses-UMCU7TSE/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_bold/', subjectcode, '_ses-UMCU7TSE_task-bairprf_MERGED_run02_bold-masked-mc-warp.nii']};

    data = {};
    for i = 1:length(nifti)
        data{i} = niftiread (nifti{i});
    end
    
    tsv = {['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-01_events.tsv'],
            ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/events/', subjectcode, '_', session, '_task-bairprf_run-02_events.tsv']};
   
    % Run1 
    events1 = tdfread(tsv{1});
    nrscans1 = size(data{1}, 4);
    onset1 = round(events1.onset / tr);
    baseline1 = find(events1.trial_type == 255);
    
    desmat1 = zeros(nrscans1, 1);   
    desmat1(onset1) = 1;
    desmat1(baseline1) = 0;
    
    desmat1 = desmat1(1:nrscans1, 1);

    % Run2
    events2 = tdfread(tsv{2});
    nrscans2 = size(data{2}, 4);
    onset2 = round(events2.onset / tr);
    baseline2 = find(events2.trial_type == 255);
    
    desmat2 = zeros(nrscans2, 1);   
    desmat2(onset2) = 1;
    desmat2(baseline2) = 0;
    
    desmat2 = desmat2(1:nrscans2, 1);
       
    % Combine design matrices
    design = {};
    design = {desmat1, desmat2};
    
    figuredir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method, '/GLM_figures'];
    
    [results_GLMdenoise, data_GLMdenoise] = GLMdenoisedata(design,data,stimdur,tr, [], [] , [], figuredir);

    disp('%%%%%%%%%%%% Done: UMCU 7T (SE) %%%%%%%%%%%%')
   
    % Save whole results structure
    output_dir = ['/Fridge/users/margriet/projects/prf/analyzeprf/results_glmdenoise/umcu/', subjectcode, '/', session, '/', method];
    cd (output_dir)
    save('data_preGLMdenoise.mat', 'data')    
    save('data_postGLMdenoise.mat', 'data_GLMdenoise')                             % Save results struct in output_dir
    save('results_GLMdenoise.mat', '-struct', 'results_GLMdenoise')            % Save results struct in output_dir
    cd ('/home/margriet/tools/prf/matlab')
    
end

%% Visualize improvements

% Voxel coordinates
xx = 77;
yy = 90;
zz = 2;

xx = 50;
yy = 32;
zz = 6;


figure(1)
data_pre = flatten(data{1}(xx, yy, zz, :));
data_post = flatten(data_GLMdenoise{1}(xx, yy, zz, :));
plot (data_pre, 'b')
hold on
plot (data_post, 'r')
hold off
xlabel ('Time (scans)')
ylabel ('BOLD signal')
legend ('Original signal', 'GLMdenoised signal')

%% End




