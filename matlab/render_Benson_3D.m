%%%%% RENDER BENSON ATLAS %%%%%

clear;
close all;
addpath(genpath('/Fridge/users/dora/github/ecogBasicCode/render/'))
addpath(genpath('/Fridge/users/margriet/projects/prf/analyzeprf/scripts/'))
addpath (genpath('/Fridge/users/giulio/github/ccep/scripts/'))
addpath (genpath('/Fridge/users/giulio/github/fieldtrip'))

%% Render Benson Areas with electrodes

dataRootPath = '/Fridge/users/dora/ccep/dataBIDS/';

subjects = {'chaam','chaam'};
hemi_cap = {'R','L'};
hemi_small = {'r','l'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','Freesurfer',['sub-' subj],'surf',...
        [hemi_small{s} 'h.benson14_varea.mgz']);
    surface_labels_B = MRIread(surface_labels_name);
    vert_label = surface_labels_B.vol(:);

    % create a colormap for the labels
    cmap = lines(max(vert_label));
    cmap(1,:) = [0.929000000000000,0.694000000000000,0.125000000000000];
    cmap(3,:) = [0.635000000000000,0.0780000000000000,0.184000000000000];
    cmap(7,:) = [0.17, 0.23, 0.82];
    cmap(9,:) = [1, 0.4, 0.7];
    cmap(10,:) = [1, 0, 0.5];
    cmap(11,:) = [0.4, 0.4, 0];
    cmap(12,:) = [0, 0.4, 0.2];
    
%     % electrode locations name:
%     dataLocName = [dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_acq-corrected_electrodes.tsv'];
%     % load electrode locations
%     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%     elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Area_Names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   

%         % make sure electrodes pop out
%         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
%         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
%         ecog_Label(els,30,12) % add electrode positions

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/BensonAreas_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end
   

