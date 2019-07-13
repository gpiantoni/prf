clear;

%%

subjectcode = 'sub-visual01'; 
subjectname = subjectcode (5:end);


%% 3TMB

% ========= 3TMB PRE ========= % 
pre3TMB = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/unfiltered/R2_masked0clust.nii']);
% ========= 3TMB POST ========= % 
post3TMB = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/final/R2_masked0clust.nii']);

nrnodes = size(pre3TMB,1) * size(pre3TMB,2) * size(pre3TMB,3);
pre3TMB = reshape (pre3TMB, [nrnodes,1]);
post3TMB = reshape (post3TMB, [nrnodes,1]);

% FIND BAD DATA IN PRE
indexInf_pre = find(isinf(pre3TMB));
indexNaN_pre = find(isnan(pre3TMB));
indexZero_pre = find (pre3TMB == 0);

% FIND BAD DATA IN POST
indexInf_post = find(isinf(post3TMB));
indexNaN_post = find(isnan(post3TMB));
indexZero_post = find (post3TMB == 0);

index_rm = cat (1, indexInf_pre, indexNaN_pre, indexZero_pre, indexInf_post, indexNaN_post, indexZero_post);

% ========= REMOVE DATA ========= % 
pre3TMB(unique(index_rm)) = [];
post3TMB(unique(index_rm)) = [];

% OUTPUT
prepost3TMB = [pre3TMB, post3TMB];

filename = ['/home/margriet/R/data/prepostFiltering/', subjectname, '_3TMB_R2_prepostFiltering'];
csvwrite (filename, prepost3TMB)

%% 7TGE

% ========= 7TGE PRE ========= % 
pre7TGE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/unfiltered/R2_masked0clust.nii']);
% ========= 7TGE POST ========= % 
post7TGE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/final/R2_masked0clust.nii']);

nrnodes = size(pre7TGE,1) * size(pre7TGE,2) * size(pre7TGE,3);
pre7TGE = reshape (pre7TGE, [nrnodes,1]);
post7TGE = reshape (post7TGE, [nrnodes,1]);

% FIND BAD DATA IN PRE
indexInf_pre = find(isinf(pre7TGE));
indexNaN_pre = find(isnan(pre7TGE));
indexZero_pre = find (pre7TGE == 0);

% FIND BAD DATA IN POST
indexInf_post = find(isinf(post7TGE));
indexNaN_post = find(isnan(post7TGE));
indexZero_post = find (post7TGE == 0);

index_rm = cat (1, indexInf_pre, indexNaN_pre, indexZero_pre, indexInf_post, indexNaN_post, indexZero_post);

% ========= REMOVE DATA ========= % 
pre7TGE(unique(index_rm)) = [];
post7TGE(unique(index_rm)) = [];

% OUTPUT
prepost7TGE = [pre7TGE, post7TGE];

filename = ['/home/margriet/R/data/prepostFiltering/', subjectname, '_7TGE_R2_prepostFiltering'];
csvwrite (filename, prepost7TGE)

%% 7TSE

% ========= 7TGE PRE ========= % 
pre7TSE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/unfiltered/R2_masked0clust.nii']);
% ========= 7TGE POST ========= % 
post7TSE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/final/R2_masked0clust.nii']);

nrnodes = size(pre7TSE,1) * size(pre7TSE,2) * size(pre7TSE,3);
pre7TSE = reshape (pre7TSE, [nrnodes,1]);
post7TSE = reshape (post7TSE, [nrnodes,1]);

% FIND BAD DATA IN PRE
indexInf_pre = find(isinf(pre7TSE));
indexNaN_pre = find(isnan(pre7TSE));
indexZero_pre = find (pre7TSE == 0);

% FIND BAD DATA IN POST
indexInf_post = find(isinf(post7TSE));
indexNaN_post = find(isnan(post7TSE));
indexZero_post = find (post7TSE == 0);

index_rm = cat (1, indexInf_pre, indexNaN_pre, indexZero_pre, indexInf_post, indexNaN_post, indexZero_post);

% ========= REMOVE DATA ========= % 
pre7TSE(unique(index_rm)) = [];
post7TSE(unique(index_rm)) = [];

% OUTPUT
prepost7TSE = [pre7TSE, post7TSE];

filename = ['/home/margriet/R/data/prepostFiltering/', subjectname, '_7TSE_R2_prepostFiltering'];
csvwrite (filename, prepost7TSE)

%% End