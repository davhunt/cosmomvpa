function [] = summarize_results()

mask_fn = '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/group_mask/group_mask_thr80pct.nii.gz';
for athresh = ["285", "357"]
  if strcmp(athresh, "285")
    tthresh = "1323";
  elseif strcmp(athresh, "357")
    tthresh = "1721";
  end
  for conn = ["1", "2", "3"]
    dirname = ['thresholds_' char(athresh) '_' char(tthresh) '_' char(conn) 'con'];
    for stat_map = ["accuracy", "tstat"]
      %for cluster_stat = ["tfce", "maxsize", "maxsum"]
      for cluster_stat = ["tfce"]
        fname = [char(stat_map) '_' char(cluster_stat) '_group_clstr_corrected_pFWE.nii.gz'];
        %p_map_nii = load_untouch_nii(['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_group_results/' dirname '/' fname]);
        p_map_ds = cosmo_fmri_dataset(['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_group_results/' dirname '/' fname],'mask',mask_fn);
        p_map = p_map_ds.samples;
        %fprintf('stat %s, clst stat %s, thr %s, conn %s: at p<0.05: %d, at p<0.01: %d, at p<0.005: %d, and at p<0.001: %d\n', stat_map, cluster_stat, athresh, conn, sum(p_map<0.05), sum(p_map<0.01), sum(p_map<0.005), sum(p_map<0.001));
        %fprintf('stat %s, clst stat %s, thr %s, conn %s: at p<0.05: %d, at p<0.01: %d, at p<0.005: %d, at p<0.001: %d, and at p<0.0005: %d\n', stat_map, cluster_stat, athresh, conn, sum(p_map<0.05), sum(p_map<0.01), sum(p_map<0.005), sum(p_map<0.001), sum(p_map<0.0005));
        fprintf('stat %s, clst stat %s, thr %s, conn %s: at p<0.05: %d, at p<0.01: %d, at p<0.005: %d, at p<0.001: %d, at p<0.0005: %d, and at p<0.0001: %d\n', stat_map, cluster_stat, athresh, conn, sum(p_map<0.05), sum(p_map<0.01), sum(p_map<0.005), sum(p_map<0.001), sum(p_map<0.0005), sum(p_map<0.0001));
        fprintf('\n');
        %fprintf('Significant voxels at p<0.005 FWE: %d, and at p<0.001 FWE: %d\n', sum(p_map<0.005), sum(p_map<0.001));
      end
    end
  end
end
