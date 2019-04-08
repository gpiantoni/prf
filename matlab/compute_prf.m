function compute_prf(nifti, output_dir, threshold)
%COMPUTE_PRF

hdr = niftiinfo(nifti{1});
n_volumes = hdr.ImageSize(3);
TR = hdr.PixelDimensions(4);

images = read_bair_stimuli(n_volumes, TR);


disp('Loading images')
for i = 1:length(nifti)
    nii{i} = niftiread(nifti{i});
end

first_nii = nii{1}(:, :,:, 1);
vxs = find(first_nii > threshold);

n_dim = size(nii{1});
n_vox = prod(n_dim(1:3));
fprintf('Selecting %d out of %d voxels (%.2f%%)\n', length(vxs), n_vox, length(vxs)/ n_vox *100)

results = analyzePRF({images, images}, nii, TR, struct('seedmode',[0 1],'display','off', 'wantglmdenoise', 0, 'vxs', vxs));


hdr.ImageSize = hdr.ImageSize(1:3);
hdr.PixelDimensions = hdr.PixelDimensions(1:3);
hdr.Datatype = 'double';

fields = {'ang', 'ecc', 'expt', 'rfsize', 'R2', 'gain'};
for i = 1:length(fields)
    niftiwrite(results.(fields{i}), fullfile(output_dir, [fields{i}  '.nii']), hdr);
end

