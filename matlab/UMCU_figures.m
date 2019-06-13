%%%%%%%%%%%% SCRIPT analyzePRF UMCU data %%%%%%%%%%%%

addpath /Fridge/users/margriet/projects/prf/analyzeprf/scripts/utilities
addpath /home/margriet/tools/prf/matlab

%% Specify subject code

subjectcode = 'sub-visual05';               % Enter subject code
session = 'ses-UMCU7TGE';

%% fMRI timeseries

% Load in data
if session == 'ses-UMCU3TMB'
    run1 = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode '_', session, '_task-bairprf_run-01_bold-rwm.nii'];
    run2 = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode '_', session, '_task-bairprf_run-02_bold-rwm.nii'];
    nifti_run1 = niftiread (run1);
    nifti_run2 = niftiread (run2);
else
    run1 = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-01_bold/', subjectcode '_', session, '_task-bairprf_run-01_bold-masked-mc-warp.nii'];
    run2 = ['/Fridge/users/margriet/subjects/bids_umcupreproc/', subjectcode, '/', session, '/', subjectcode, '_', session, '_task-bairprf_run-02_bold/', subjectcode '_', session, '_task-bairprf_run-02_bold-masked-mc-warp.nii'];
    nifti_run1 = niftiread (run1);
    nifti_run2 = niftiread (run2);
end

% Reshape nifti to 2D matrix
dim = size (nifti_run1);       % [fx X fy X fz X dynamics]
fx = dim(1);
fy = dim(2);
fz = dim(3);
nrnodes = fx*fy*fz;          % To convert 4D to 2D = 114688

nrscans_run1 = size(nifti_run1, 4);       % dynamics = 224
nrscans_run2 = size(nifti_run2, 4);

nifti_run1 = reshape(nifti_run1, [nrnodes, nrscans_run1]);  % DIM = [114688 X 224]
nifti_run2 = reshape(nifti_run2, [nrnodes, nrscans_run1]);

data_UMCU = {nifti_run1, nifti_run2};

% FIGURE: Plot times series per voxel
voxelnr = 2200;                % Select voxel
vxstr = num2str (voxelnr);

temp = cellfun(@(x) squeeze(x(voxelnr,:)),data_UMCU,'UniformOutput',0);
figure (1)
hold on;
set(gcf,'Units','points','Position',[100 100 600 150]);
plot(cat(2,temp{:}),'r-');
straightline(224*(1:2)+.5,'v','g-');
hold off
xlabel('TR');
ylabel('BOLD signal');
ax = axis;
axis([.5 600+.5 ax(3:4)]);
title(['Time-series data: voxel ', vxstr]);

%% Inspect the stimuli

str_data_UMCU = {run1, run2};

images = {};
% 1st run
run = 1;
hdr = niftiinfo(str_data_UMCU{1});                % read nifti header
TR = hdr.PixelDimensions(4);                      % 850 ms
n_volumes = hdr.ImageSize(4);                     % no. of dynamics
images_run1 = read_bair_stimuli(subjectcode, session, run, n_volumes, TR);

% 2nd run
run = 2;
hdr = niftiinfo(str_data_UMCU{2});                % read nifti header
TR = hdr.PixelDimensions(4);                      % 850 ms
n_volumes = hdr.ImageSize(4);                     % no. of dynamics
images_run2 = read_bair_stimuli(subjectcode, session, run, n_volumes, TR);

images = {images_run1, images_run2};

% FIGURE: Plot stimuli for 3 time points
figure(2)
set(gcf,'Units','points','Position',[100 100 700 300]);
for p=1:3
  subplot(1,3,p); 
  hold on;
  num = 30+2*p;
  imagesc(images{1}(:,:,num),[0 1]);
  axis image tight;
  set(gca,'YDir','reverse');
  colormap(gray);
  title(sprintf('Image number %d',num));
end
hold off

%% Visualize pRF results

visualangle3T = 13.0350;
visualangle7T = 12.8560;

if session == 'ses-UMCU3TMB'
    cfactor = visualangle3T/100;
else
    cfactor = visualangle7T/100;
end

% Load in results
ang = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/ang.nii']);
ecc = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/ecc.nii']);
rfsize = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/rfsize.nii']);

% Reshape results
ang = reshape(ang, [nrnodes, 1]);
ecc = reshape(ecc, [nrnodes, 1]);
rfsize = reshape(rfsize, [nrnodes, 1]);

voxmin = 53580;
voxmax = 53585;
strvxmin = num2str(voxmin);
strvxmax = num2str(voxmax);

% Results voxels 53593 - 53598
results_UMCU.ang = ang(voxmin:voxmax);
results_UMCU.ecc = ecc(voxmin:voxmax);
results_UMCU.rfsize = rfsize(voxmin:voxmax);

% FIGURE: Results
figure (3); 
hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
title (['pRF parameters for voxels: ', strvxmin, '-', strvxmax])
drawrectangle(0,0,10,10,'k-'); 
cmap = jet(size(results_UMCU.ang,1));
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-10:2:10,'YTick',-10:2:10);
for  p=1:size(results_UMCU.ang,1)
  xpos = squeeze(results_UMCU.ecc(p) * cos(results_UMCU.ang(p)/180*pi) * cfactor);
  ypos = squeeze(results_UMCU.ecc(p) * sin(results_UMCU.ang(p)/180*pi) * cfactor);
  angle = squeeze(results_UMCU.ang(p)/180*pi);
  sd = squeeze(results_UMCU.rfsize(p) * cfactor);
  h = drawellipse(xpos,ypos,angle,2*sd,2*sd);  % circle at +/- 2 pRF sizes
  set(h,'Color',cmap(p,:),'LineWidth',2);
  set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
end
hold off

%% Model fit

% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = load('hrf');                  % HRF that was used in the model
hrf = hrf.hrf;
degs = [2,2];  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(images)
  stimulusPP{p} = squish(images{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

polymatrix = {};
for p=1:length(degs)
  polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data_UMCU{p},2),0:degs(p)));
end

% Saved results
results = load (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/', session, '/separate_bairprf/results.mat']);
results_UMCU.params = results.params;

R2 = results.R2;
R2 = reshape(R2, [nrnodes, 1]);

% Only include voxels above threshold
index = find (all(~isnan(R2),2));

data_UMCU{1} = data_UMCU{1}(index,:);
data_UMCU{2} = data_UMCU{2}(index,:);

R2 = R2(index);

% Find max R2
[maxR2, indexR2] = max(R2);

% Which voxel should we inspect?  
vx = indexR2;

datats = {};
modelts = {};
for p=1:length(data_UMCU)
  datats{p} =  polymatrix{p}*data_UMCU{p}(vx,:)';
  modelts{p} = polymatrix{p}*modelfun(results_UMCU.params(1,:,vx),stimulusPP{p});
end

% FIGURE: Modelfit
figure (4); 
hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'r-');
plot(cat(1,modelts{:}),'b-');
straightline(nrscans_run1*(1:2)+.5,'v','g-');
xlabel('Time (s)');
ylabel('BOLD signal');
ax = axis;
axis([.5 500+.5 ax(3:4)]);
title('Time-series data - Modelfit');
hold off

%%
