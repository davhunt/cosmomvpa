function [] = cosmo_groupstats(subfolder, outdir, nIter, opt)

arguments
    subfolder
    outdir
    nIter
    opt.output_string = ''
    opt.cluster_stat = 'tfce'
    opt.connectivity = 1
    opt.null_data_mat = 'none'
    opt.nNulls_per_subject = 500
    opt.dh = 0.1
    opt.cluster_threshold = 2.85
    opt.stat_map = 'accuracy'
end

% nIter: number of permuted group average iterations used to generate the null distribution
% subfolder: folder under Complex_seq_analysis/XXX/MVPA to look for perm and results folders
% outdir: folder (full path) to save results to
% output_string: optional string to prepend to output file names (default: '')
% cluster_stat: 'tfce', 'maxsize', or 'maxsum'; cluster statistic to generate p-value map with (default: 'tfce')
% connectivity: how to define contiguous voxels, 1 for just faces, 2 for face or edge, 3 for face edge or vertex (default: 1)
% null_data_mat: path to null data .mat file, or 'none' to generate from "perm" folder (default: none)
% nNulls_per_subject: how many nulls from each subject are available to randomly choose from to build group map (default: 500)
% dh: size of integration step dh (default: 0.1)
% cluster_threshold: cluster-forming threshold for "maxsize" or "maxsum" argument to "cosmo_measure_clusters" (default: 2.85) 
%                        keep in mind if using t-statistic threshold should be a t-value (i.e. 1.721 for p < 0.05 and df = 21, or 1.323 for p < 0.1 which is reasonable)
% stat_map: 'accuracy' or 'tstat' (default: 'accuracy')

%% sample call:
%% cosmo_groupstats(subfolder, outdir, nIter, cluster_stat='tfce', connectivity=1, null_data_mat='/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_group_results/null_data.mat', nNulls_per_subject=500, dh=0.1, cluster_threshold=10.0, stat_map='accuracy')

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
subj_list = [1 2 3 4 5 6 7 8 9 10 13 14 15 16 17 18 19 20 21 22 23 24];
% Could modify this
% subj_list = [2 4 14 17 18 19];
nSubjects = numel(subj_list);
output_string = opt.output_string;
nNulls_per_subject = opt.nNulls_per_subject;
null_data_mat = opt.null_data_mat;
connectivity = opt.connectivity;
cluster_stat = opt.cluster_stat;
dh = opt.dh;
cluster_threshold = opt.cluster_threshold;
stat_map = opt.stat_map;
%outdir = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_group_results/better_mask';
% ---------------------------

%mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask_thr_85pct_bin_bps.nii.gz';
mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask/group_mask_thr80pct.nii.gz';
%mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask_thr_85_new.nii.gz';
% ^ A group mask consisting of all voxels present in >= 85% of the individual subjects'' masks
% at least 19 / 22 subjects

assert(isfile(mask_fn), 'can''t find mask file');
assert(ismember(cluster_stat, {'tfce','maxsize','maxsum'}), 'cluster_stat must be tfce or maxsize');
assert(ismember(stat_map, {'accuracy','tstat'}), 'stat_map must be accuracy or tstat');

subjs = {};
for i = subj_list
    sub = sprintf('%02d', i);
    %subjs{end+1} = ['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/TDT_results/TDT_results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz'];
    subjs{end+1} = ['/N/project/brainevo/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/TDT_results/TDT_results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz'];
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

%ds_group = cosmo_stack(ds_all);


if strcmp(null_data_mat, 'none') == 1
    null_data = cell(nSubjects,1);
    for s=1:nSubjects
        fprintf('Loading nulls for subject %d/%d\n',s, numel(subj_list));
        fnames = cell(nNulls_per_subject,1);
        for i = 1:nNulls_per_subject
          %fnames{i} = sprintf('/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/0%02d/MVPA/%s/perm/perm%03d/results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz', subj_list(s), subfolder, i);
          fnames{i} = sprintf('/N/project/brainevo/Complex_seq_analysis/0%02d/MVPA/%s/perm/perm%03d/results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz', subj_list(s), subfolder, i);
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
nbrhood = cosmo_cluster_neighborhood(template_ds, 'fmri', connectivity); % connectivity 1 2 or 3
P = numel(nbrhood.neighbors);
max_neighbors = max(cellfun(@numel, nbrhood.neighbors));
nbrhood_mat = zeros(max_neighbors,P,'uint32');
for p = 1:P
  idxs = nbrhood.neighbors{p}(:)';
  m = numel(idxs);
  nbrhood_mat(1:m, p) = idxs;
end

%% Cluster statistic map from observed group mean accuracy map
if strcmp(stat_map, 'accuracy')
    group_obs = mean(obs_data,1); % 1 x nVoxels
elseif strcmp(stat_map, 'tstat')
    %obs_data = obs_data ./ 100; % normalize to [-.5, .5] for t-test %actually not necessary
    t_obs = onesample_t(obs_data);
    group_obs = t_obs;
end

%tfce_obs  = cosmo_tfce(group_obs, template_ds.fa, template_ds.a); % function doesn't exist?
if strcmp(cluster_stat, 'tfce')
    clusterstat_obs = cosmo_measure_clusters(group_obs, nbrhood_mat, cluster_stat, 'dh', dh);
elseif strcmp(cluster_stat, 'maxsize') || strcmp(cluster_stat, 'maxsum')
    clusterstat_obs = cosmo_measure_clusters(group_obs, nbrhood_mat, cluster_stat, 'threshold', cluster_threshold);
end
%tfce_null = cosmo_measure_clusters(group_null, nbrhood_mat, 'tfce', 'dh', dh);

%% Build null distribution of group mean accuracy maps
null_max_clusterstat = zeros(nIter,1);

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
    if strcmp(stat_map, 'accuracy')
        group_null = mean(null_sample,1); % 1 x nVoxels
    elseif strcmp(stat_map, 'tstat')
        %null_sample = null_sample ./ 100; % normalize to [-.5, .5] for t-test %actually not necessary
        t_null = onesample_t(null_sample);
        group_null = t_null;
    end
    
    % TFCE, or maxsize
    if strcmp(cluster_stat, 'tfce')
        clusterstat_null = cosmo_measure_clusters(group_null, nbrhood_mat, 'tfce', 'dh', dh); 
    elseif strcmp(cluster_stat, 'maxsize') || strcmp(cluster_stat, 'maxsum')
        clusterstat_null = cosmo_measure_clusters(group_null, nbrhood_mat, cluster_stat, 'threshold', cluster_threshold);
    end
    
    % store maximum TFCE, maxsize, or maxsum value across voxels to build null distribution of max vals
    null_max_clusterstat(it) = max(clusterstat_null);
end

%% Compute voxelwise FWE-corrected p-values
p_map = zeros(1,nVoxels);
for v = 1:nVoxels
    p_map(v) = (sum(null_max_clusterstat >= clusterstat_obs(v)) + 1) / (nIter + 1);
end

%% Save results

if ~isfolder(outdir)
  mkdir(outdir)
end

ds_out = template_ds;
ds_out.samples = clusterstat_obs;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_observed.nii.gz']);

ds_out.samples = p_map;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE.nii.gz']);

ds_out.samples = p_map < 0.05;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE05_clusters.nii.gz']);

ds_out.samples = p_map < 0.01;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE01_clusters.nii.gz']);

fprintf('Done. Significant voxels at p<0.05 FWE: %d, and at p<0.01 FWE: %d\n', sum(p_map<0.05), sum(p_map<0.01));
fprintf('Significant voxels at p<0.005 FWE: %d, and at p<0.001 FWE: %d\n', sum(p_map<0.005), sum(p_map<0.001));

end

function tvec = onesample_t(X)
    % Compute t-map on X
    % X: [nSubjects x nVoxels]
    % tvec: [1 x nVoxels]
    n  = size(X,1);
    mu = mean(X,1);
    sd = std(X,0,1);           % sample std (dof = n-1)
    se = sd ./ sqrt(n);
    tvec = zeros(1,size(X,2));
    nz = se>0 & isfinite(se);
    tvec(nz) = mu(nz) ./ se(nz);
    % if se==0 leave t=0
end
