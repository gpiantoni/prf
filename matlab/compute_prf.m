function compute_prf(nifti, output_dir, threshold)

addpath(genpath('/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF'))

threshold = 2000;
nifti =  {'/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-01_preproc.nii.gz',
    '/Fridge/users/giovanni/projects/margriet/bids_nyu/derivatives/preprocessed/sub-visual01/ses-UMCU3TMB/sub-visual01_ses-UMCU3TMB_task-bairprf_run-02_preproc.nii.gz'};
output_dir = '/Fridge/users/margriet/projects/analysis_code/code_Kay_analyzePRF/results_analyzePRF/01';

images = read_bair_stimuli();

hdr = niftiinfo(nifti{1});
TR = hdr.PixelDimensions(4);

disp('Loading images')
for i = 1:length(nifti)
    nii{i} = niftiread(nifti{i});
end

nii_cat = cat(4, nii{:});
mean_nii = mean(nii_cat, 4);
vxs = find(mean_nii > threshold);

n_dim = size(nii{1});
n_vox = prod(n_dim(1:3));
fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)

results = analyzePRF({images, images}, nii, TR, struct('seedmode',[0 1],'display','off', 'wantglmdenoise', 0, 'vxs', find(vxs)));


hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';

fields = {'ang', 'ecc', 'expt', 'rfsize', 'R2', 'gain'};
for i = 1:length(fields)
    niftiwrite(results.(fields{i}), fullfile(output_dir, [fields{i}  '.nii']), hdr);
end

