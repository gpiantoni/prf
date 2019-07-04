% % =========  SCRIPT pRF R2 values per session  ========= % 
%
%%
clear;

%%

subjectcode = 'sub-visual01'; 
method = 'final';

% ========= 3TMB ========= % 
R2_3TMB = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU3TMB/', method, '/R2_masked0clust.nii']);

nrnodes = size(R2_3TMB,1) * size(R2_3TMB,2) * size(R2_3TMB,3);

R2_3TMB = reshape (R2_3TMB, [nrnodes,1]);
indexInf = find(isinf(R2_3TMB));
indexNaN = find(isnan(R2_3TMB));
indexZero = find (R2_3TMB == 0);
index_rm = cat (1, indexInf, indexNaN, indexZero);
R2_3TMB(unique(index_rm)) = [];

R2_total = size(R2_3TMB(R2_3TMB>0),1);
R2_5 = size(R2_3TMB(R2_3TMB>5),1);
R2_10 = size(R2_3TMB(R2_3TMB>10),1);
R2_15 = size(R2_3TMB(R2_3TMB>15),1);
R2_20 = size(R2_3TMB(R2_3TMB>20),1);
R2_25 = size(R2_3TMB(R2_3TMB>25),1);
R2_40 = size(R2_3TMB(R2_3TMB>40),1);

val_3TMB = [R2_total, R2_5, R2_10, R2_15, R2_20, R2_25, R2_40];
per_3TMB = [(R2_total/R2_total)*100, (R2_5/R2_total)*100, (R2_10/R2_total)*100, (R2_15/R2_total)*100, (R2_20/R2_total)*100, (R2_25/R2_total)*100, (R2_40/R2_total)*100];

% ========= 7TGE ========= % 
R2_7TGE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TGE/', method, '/R2_masked0clust.nii']);

nrnodes = size(R2_7TGE,1) * size(R2_7TGE,2) * size(R2_7TGE,3);
R2_7TGE = reshape (R2_7TGE, [nrnodes,1]);

indexInf = find(isinf(R2_7TGE));
indexNaN = find(isnan(R2_7TGE));
indexZero = find (R2_7TGE == 0);
index_rm = cat (1, indexInf, indexNaN, indexZero);
R2_7TGE(unique(index_rm)) = [];

R2_total = size(R2_7TGE(R2_7TGE>0),1);
R2_5 = size(R2_7TGE(R2_7TGE>5),1);
R2_10 = size(R2_7TGE(R2_7TGE>10),1);
R2_15 = size(R2_7TGE(R2_7TGE>15),1);
R2_20 = size(R2_7TGE(R2_7TGE>20),1);
R2_25 = size(R2_7TGE(R2_7TGE>25),1);
R2_40 = size(R2_7TGE(R2_7TGE>40),1);

val_7TGE = [R2_total, R2_5, R2_10, R2_15, R2_20, R2_25, R2_40];
per_7TGE = [(R2_total/R2_total)*100, (R2_5/R2_total)*100, (R2_10/R2_total)*100, (R2_15/R2_total)*100, (R2_20/R2_total)*100, (R2_25/R2_total)*100, (R2_40/R2_total)*100];

% ========= 7TSE ========= % 
R2_7TSE = niftiread (['/Fridge/users/margriet/projects/prf/analyzeprf/results/umcu/', subjectcode, '/ses-UMCU7TSE/', method, '/R2_masked0clust.nii']);

nrnodes = size(R2_7TSE,1) * size(R2_7TSE,2) * size(R2_7TSE,3);
R2_7TSE = reshape (R2_7TSE, [nrnodes,1]);

indexInf = find(isinf(R2_7TSE));
indexNaN = find(isnan(R2_7TSE));
indexZero = find (R2_7TSE == 0);
index_rm = cat (1, indexInf, indexNaN, indexZero);
R2_7TSE(unique(index_rm)) = [];

R2_total = size(R2_7TSE(R2_7TSE>0),1);
R2_5 = size(R2_7TSE(R2_7TSE>5),1);
R2_10 = size(R2_7TSE(R2_7TSE>10),1);
R2_15 = size(R2_7TSE(R2_7TSE>15),1);
R2_20 = size(R2_7TSE(R2_7TSE>20),1);
R2_25 = size(R2_7TSE(R2_7TSE>25),1);
R2_40 = size(R2_7TSE(R2_7TSE>40),1);

val_7TSE = [R2_total, R2_5, R2_10, R2_15, R2_20, R2_25, R2_40];
per_7TSE = [(R2_total/R2_total)*100, (R2_5/R2_total)*100, (R2_10/R2_total)*100, (R2_15/R2_total)*100, (R2_20/R2_total)*100, (R2_25/R2_total)*100, (R2_40/R2_total)*100];

%% Summary per subject
 
val = [val_3TMB; val_7TGE; val_7TSE];
per = [per_3TMB; per_7TGE; per_7TSE];


