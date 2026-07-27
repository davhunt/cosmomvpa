function [] = cosmo_sigclusters(subfolder, outdir, subj, opt)

arguments
    subfolder
    outdir
    subj
    opt.output_string = 'tstat_tfce_'
    opt.cluster_stat = 'tfce'
    opt.connectivity = 1
    opt.null_data_mat = 'none'
    opt.nNulls_per_subject = 1000
    opt.dh = 0.1
    opt.cluster_threshold = 2.831
end

% Script to perform cluster correction and generate significant clusters maps at the per-subject level from a subject's true-labelled decoding accuracy maps (compare to cosmo_groupstats.m)
%
% subfolder: folder under Complex_seq_analysis/XXX/MVPA to look for perm and results folders ('breakpoint_rts_unsmoothed_native_antsmask_new_bps')
% outdir: folder under Complex_seq_analysis/XXX/MVPA/<subfolder> to save results to
% subj: subject number
% cluster_stat: 'tfce', 'maxsize', or 'maxsum'; cluster statistic to generate p-value map with (default: 'tfce')
% connectivity: how to define contiguous voxels, 1 for just faces, 2 for face or edge, 3 for face edge or vertex (default: 1)
% null_data_mat: path to null data .mat file, or 'none' to generate from "perm" folder (default: none)
% nNulls_per_subject: how many nulls from each subject are available to randomly choose from to build group map (default: 500)
% dh: size of integration step dh (default: 0.1)
% cluster_threshold: cluster-forming threshold for "maxsize" or "maxsum" argument to "cosmo_measure_clusters" (default: 2.85) 
%                        keep in mind if using t-statistic threshold should be a t-value (i.e. 1.721 for p < 0.05 and df = 21, or 1.323 for p < 0.1 which is reasonable)




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
output_string = opt.output_string;
nNulls_per_subject = opt.nNulls_per_subject;
null_data_mat = opt.null_data_mat;
connectivity = opt.connectivity;
cluster_stat = opt.cluster_stat;
dh = opt.dh;
cluster_threshold = opt.cluster_threshold;
%outdir = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_group_results/better_mask';
% ---------------------------

% get mask-- might as well stay in MNI space?
if exist('/N/slate/brainevo/Implicit_Learning/CoSMoMVPA/mvpa') == 7 % slate
    mask_fn = ['/N/project/brainevo/Complex_seq_analysis/0' num2str(subj, '%.2d') '/MVPA/breakpoint_rts_unsmoothed_native_antsmask_new_bps/mask_mni.nii'];
else
    mask_fn = ['/Volumes/data1/Complex_seq_analysis/0' num2str(subj, '%.2d') '/MVPA/breakpoint_rts_unsmoothed_native_antsmask_new_bps/mask_mni.nii'];
end


assert(isfile(mask_fn), 'can''t find mask file');
assert(ismember(cluster_stat, {'tfce','maxsize','maxsum'}), 'cluster_stat must be tfce or maxsize');

%subj = num2str(subj, '%.2d');
acc_map = ['/N/project/brainevo/Complex_seq_analysis/0' num2str(subj, '%.2d') '/MVPA/' subfolder '/TDT_results/TDT_results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz'];


ds_sub = cosmo_fmri_dataset(acc_map,'mask',mask_fn);
obs_data = ds_sub.samples; % 1 x nVoxels matrix of decoding accuracy values

%obs_data = [];
%ds_all = cell(numel(subjs),1);
%for s = 1:numel(subjs)
%    ds_all{s} = cosmo_fmri_dataset(subjs{s},'mask',mask_fn);
%    ds_all{s}.sa.targets = 1;
%    ds_all{s}.sa.chunks = s;
%    obs_data(s,:) = ds_all{s}.samples;
%end
nVoxels = size(ds_sub.samples,2);

%ds_group = cosmo_stack(ds_all);

if strcmp(null_data_mat, 'none') == 1
    null_data = zeros(nNulls_per_subject, nVoxels);
    for i = 1:nNulls_per_subject
        fname = sprintf('/N/project/brainevo/Complex_seq_analysis/0%02d/MVPA/%s/perm/perm%03d/results_avg_allreps/res_accuracy_minus_chance_mni.nii.gz', subj, subfolder, i);
        ds_tmp = cosmo_fmri_dataset(fname,'mask',mask_fn);
        null_data(i,:) = ds_tmp.samples;
    end
else
    load(null_data_mat); % save a few minutes loading null_data from a pre-saved .mat file
end


rng(0,'twister'); % reproducibility
%% Construct nbrhood_mat matrix from nbrhood to feed into cosmo_measure_clusters 
nbrhood = cosmo_cluster_neighborhood(ds_sub, 'fmri', connectivity); % connectivity 1 2 or 3
P = numel(nbrhood.neighbors);
max_neighbors = max(cellfun(@numel, nbrhood.neighbors));
nbrhood_mat = zeros(max_neighbors,P,'uint32');
for p = 1:P
  idxs = nbrhood.neighbors{p}(:)';
  m = numel(idxs);
  nbrhood_mat(1:m, p) = idxs;
end

if strcmp(cluster_stat, 'tfce')
    clusterstat_obs = cosmo_measure_clusters(obs_data, nbrhood_mat, cluster_stat, 'dh', dh);
elseif strcmp(cluster_stat, 'maxsize') || strcmp(cluster_stat, 'maxsum')
    clusterstat_obs = cosmo_measure_clusters(obs_data, nbrhood_mat, cluster_stat, 'threshold', cluster_threshold);
end

null_max_clusterstat = zeros(nNulls_per_subject,1);
for i = 1:nNulls_per_subject
    null_map = null_data(i,:);
    if strcmp(cluster_stat, 'tfce')
        clusterstat_null = cosmo_measure_clusters(null_map, nbrhood_mat, 'tfce', 'dh', dh); 
    elseif strcmp(cluster_stat, 'maxsize') || strcmp(cluster_stat, 'maxsum')
        clusterstat_null = cosmo_measure_clusters(null_map, nbrhood_mat, cluster_stat, 'threshold', cluster_threshold);
    end
    % store maximum TFCE, maxsize, or maxsum value across voxels to build null distribution of max vals
    null_max_clusterstat(i) = max(clusterstat_null);
end


%% Compute voxelwise FWE-corrected p-values
p_map = zeros(1,nVoxels);
for v = 1:nVoxels
    p_map(v) = (sum(null_max_clusterstat >= clusterstat_obs(v)) + 1) / (nNulls_per_subject + 1);
end

%% Save results

%if ~isfolder(outdir)
%  mkdir(outdir)
%end
outdir = ['/N/project/brainevo/Complex_seq_analysis/0' num2str(subj, '%.2d') '/MVPA/' subfolder '/' outdir];
if ~isfolder(outdir)
  mkdir(outdir);
end

ds_out = ds_sub;
ds_out.samples = clusterstat_obs;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_observed.nii.gz']);

ds_out.samples = p_map;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE.nii.gz']);

ds_out.samples = p_map < 0.05;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE05_clusters.nii.gz']);

ds_out.samples = p_map < 0.01;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE01_clusters.nii.gz']);

ds_out.samples = p_map < 0.005;
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE005_clusters.nii.gz']);

ds_out.samples = p_map < 0.001; % This is the smallest possible p value for nNulls_per_subject = 1000 (1 / 1001)
cosmo_map2fmri(ds_out, [ outdir '/' char(output_string) 'group_clstr_corrected_pFWE001_clusters.nii.gz']);


fprintf('Subject %d done. Significant voxels at p<0.05 FWE: %d, and at p<0.01 FWE: %d\n', subj, sum(p_map<0.05), sum(p_map<0.01));
fprintf('Significant voxels at p<0.005 FWE: %d, and at p<0.001 FWE: %d\n', sum(p_map<0.005), sum(p_map<0.001));

end
