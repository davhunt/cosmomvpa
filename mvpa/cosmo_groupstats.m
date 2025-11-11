function [] = cosmo_groupstats(subfolder, outdir, nIter, cluster_stat, connectivity, null_data_mat)

% nIter: number of permuted group average iterations used to generate the null distribution
% subfolder: folder under Complex_seq_analysis/XXX/MVPA to look for perm and results folders
% cluster_stat: 'tfce' or 'maxsize'; cluster statistic to generate p-value map with
% connectivity: how to define contiguous voxels, 1 for just faces, 2 for face or edge, 3 for face edge or vertex
% null_data_mat: path to null data .mat file, or 'none' to generate from "perm" folder

if exist('/N/slate/brainevo/Implicit_Learning/CoSMoMVPA/mvpa') == 7 % slate
    addpath(genpath('/N/slate/brainevo/Implicit_Learning/CoSMoMVPA/mvpa'));
    addpath(genpath('/N/slate/brainevo/Implicit_Learning/matlab_nifti_tools'));
else
    addpath(genpath('/Users/lab/Downloads/CoSMoMVPA/mvpa')); % habilis
    addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools'));
end

% ---------------------------
% Parameters (edit as needed)
% ---------------------------
%nIter = 10000;          % number of group-level null iterations
subj_list = [1 2 3 4 5 6 7 8 9 10 13 14 15 16 17 18 19 20 21 22 23 24];
% Could modify this
% subj_list = [2 4 14 17 18 19];
nSubjects = numel(subj_list);
nNulls_per_subject = 500;
%outdir = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_group_results/better_mask';
% ---------------------------

mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask_thr_85pct_bin_bps.nii.gz';
%mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask_thr_85_new.nii.gz';
% ^ A group mask consisting of all voxels present in >= 85% of the individual subjects' masks
% at least 19 / 22 subjects

subjs = {};
for i = subj_list
    sub = sprintf('%02d', i);
    subjs{end+1} = ['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/TDT_results/TDT_results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz'];
end

template_ds = cosmo_fmri_dataset(subjs{1},'mask',mask_fn);

obs_data = [];
ds_all = cell(numel(subjs),1);
for s = 1:numel(subjs)
    ds_all{s} = cosmo_fmri_dataset(subjs{s},'mask',mask_fn);
    ds_all{s}.sa.targets = 1;
    ds_all{s}.sa.chunks = s;
    obs_data(s,:) = ds_all{s}.samples;
end
nVoxels = size(obs_data,2);

ds_group = cosmo_stack(ds_all);


if strcmp(null_data_mat, 'none') == 1
    for s=1:nSubjects
        fprintf('Loading nulls for subject %d/%d\n',s, numel(subj_list));
        fnames = cell(nNulls_per_subject,1);
        for i = 1:nNulls_per_subject
          fnames{i} = sprintf('/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/0%02d/MVPA/breakpoint_rts_unsmoothed_native/perm/perm%03d/results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz', subj_list(s), i);
        end
        assert(numel(fnames) == nNulls_per_subject, 'Mismatch for sub %d',subj_list(s));
        tmp = zeros(nNulls_per_subject, nVoxels);
        for i = 1:nNulls_per_subject
            ds_tmp = cosmo_fmri_dataset(fnames{i},'mask',mask_fn);
            tmp(i,:) = ds_tmp.samples;
        end
    null_data{s} = tmp;
    end
else
    load(null_data_mat); % save a few minutes loading null_data from a pre-saved .mat file
end




% Inputs:
%  obs_data   : [nSubjects x nVoxels] observed accuracy maps
%  null_data  : {nSubjects x 1}, each cell [nNulls_per_subject x nVoxels] = null maps for that subject
%  template_ds: CoSMo dataset with .fa and .a (for cosmo_tfce and saving output)
%  nIter      : number of null iterations to sample (e.g. 5000)

%nSubjects             = size(obs_data,1);
%nVoxels               = size(obs_data,2);
%nNulls_per_subject    = size(null_data{1},1);

rng(0,'twister'); % reproducibility
%% Construct nbrhood_mat matrix from nbrhood to feed into cosmo_measure_clusters 
%nbrhood = cosmo_cluster_neighborhood(template_ds);
nbrhood = cosmo_cluster_neighborhood(template_ds, 'fmri', connectivity); % connectivity 1 2 or 3
%nbrhood = cosmo_cluster_neighborhood(template_ds, 'fmri', 1); % connectivity 1 2 or 3
P = numel(nbrhood.neighbors);
max_neighbors = max(cellfun(@numel, nbrhood.neighbors));
nbrhood_mat = zeros(max_neighbors,P,'uint32');
for p = 1:P
  idxs = nbrhood.neighbors{p}(:)';
  m = numel(idxs);
  nbrhood_mat(1:m, p) = idxs;
end





%% --- Step 1. Observed group map (mean across subjects) ---
group_obs = mean(obs_data,1); % 1 x nVoxels
%tfce_obs  = cosmo_tfce(group_obs, template_ds.fa, template_ds.a);
tfce_obs = cosmo_measure_clusters(group_obs, nbrhood_mat, 'tfce', 'dh', 0.1);
%tfce_null = cosmo_measure_clusters(group_null, nbrhood_mat, 'tfce', 'dh', 0.1);

%% --- Step 2. Null distribution ---
null_max_tfce = zeros(nIter,1);

for it = 1:nIter
    if mod(it,500)==0
        fprintf('Iteration %d/%d\n', it, nIter);
    end
    
    % sample one null map per subject
    null_sample = zeros(nSubjects, nVoxels);
    for s = 1:nSubjects
        idx = randi(nNulls_per_subject); % random null index
        null_sample(s,:) = null_data{s}(idx,:);
    end
    
    % group average
    group_null = mean(null_sample,1);
    
    % TFCE
    %tfce_null = cosmo_tfce(group_null, template_ds.fa, template_ds.a);
    tfce_null = cosmo_measure_clusters(group_null, nbrhood_mat, 'tfce', 'dh', 0.1); 
    
    % store maximum TFCE value across voxels
    null_max_tfce(it) = max(tfce_null);
end

%% --- Step 3. Compute voxelwise FWE-corrected p-values ---
p_map = zeros(1,nVoxels);
for v = 1:nVoxels
    p_map(v) = (sum(null_max_tfce >= tfce_obs(v)) + 1) / (nIter + 1);
end

%% --- Step 4. Save results ---

if ~isfolder(outdir)
  mkdir(outdir)
end

ds_out = template_ds;
ds_out.samples = tfce_obs;
cosmo_map2fmri(ds_out, [ outdir '/group_TFCE_observed.nii.gz']);

ds_out.samples = p_map;
cosmo_map2fmri(ds_out, [ outdir '/group_TFCE_pFWE.nii.gz']);

ds_out.samples = p_map < 0.05;
cosmo_map2fmri(ds_out, [ outdir '/group_TFCE_pFWE05_mask.nii.gz']);

fprintf('Done. Significant voxels at p<0.05 FWE: %d\n', sum(p_map<0.05));

