function [] = generate_results()

null_data_mat='/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_group_results/null_data.mat';
subfolder = 'breakpoint_rts_unsmoothed_native_antsmask';
nIter = 10000;

%for stat_map = ["accuracy", "tstat"]
for stat_map = ["tstat", "accuracy"]
    for cluster_stat = ["tfce", "maxsize", "maxsum"]
        ostr = strcat(stat_map, "_", cluster_stat, "_");
        if strcmp(stat_map, "accuracy")
            cosmo_groupstats(subfolder, outdir, nIter, output_string=ostr, cluster_stat=cluster_stat, ...
              connectivity=1, null_data_mat=null_data_mat, nNulls_per_subject=500, dh=0.1, ...
              cluster_threshold=10.0, stat_map=stat_map)
        elseif strcmp(stat_map, "tstat")
            cosmo_groupstats(subfolder, outdir, nIter, output_string=ostr, cluster_stat=cluster_stat, ...
              connectivity=1, null_data_mat=null_data_mat, nNulls_per_subject=500, dh=0.1, ...
              cluster_threshold=1.323, stat_map=stat_map)
        end
    end
end
