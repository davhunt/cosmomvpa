function [] = generate_results()

for athresh = ["285", "357"]
  if strcmp(athresh, "285")
    ttresh = "1323";
  elseif strcmp(athresh, "357")
    ttresh = "1721";
  end
  for conn = ["1", "2", "3"]
    dirname = ['thresholds_' char(athresh) '_' char(tthresh) '_' char(conn) 'con'];
    for stat_map = ["accuracy", "tstat"]
      %for cluster_stat = ["tfce", "maxsize", "maxsum"]
      for cluster_stat = ["tfce"]
        fname = [char(stat_map) '_' char(cluster_stat) '_group_clstr_corrected_pFWE.nii.gz'];
        p_map_nii = load_untouch_nii(['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_group_results/' dirname '/' fname]);
        p_map = p_map_nii.img;
        fprintf('stat %s, clst stat %s, thr %s, conn %s: at p<0.05: %d, at p<0.01: %d, at p<0.005: %d, and at p<0.001: %d\n', stat_map, cluster_stat, athresh, conn, sum(p_map<0.05), sum(p_map<0.01), sum(p_map<0.005), sum(p_map<0.001));
        %fprintf('Significant voxels at p<0.005 FWE: %d, and at p<0.001 FWE: %d\n', sum(p_map<0.005), sum(p_map<0.001));
      end
    end
  end
end
