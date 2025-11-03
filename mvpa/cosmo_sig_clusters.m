function [] = cosmo_sig_clusters(subj, dh, subfolder, nperms)

% subj: subject number
% dh: threshold step size (For typical use cases, a value so that 100*dh is in the same order of magnitude as the range)
% subfolder: Folder under MVPA to find "perm" folder ("durations_unsmoothed_native")
% number of permuations in "perm" folder

if exist('/N/slate/brainevo/Implicit_Learning') == 7
    base_folder = ['/N/slate/brainevo/Implicit_Learning']; % If on slate
    addpath(genpath([base_folder '/matlab_nifti_tools']));
    addpath(genpath([base_folder '/CoSMoMVPA/mvpa']));
elseif exist('/Volumes/data1') == 7
    base_folder = ['/Volumes/data1']; % Or on habilis
    addpath(genpath('/Users/lab/Downloads/CoSMoMVPA/mvpa'));
    addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools'));
else
    error('Can''t find base folder, are we on habilis or the HPC?');
end

true_data_nii = cosmo_fmri_dataset([base_folder '/Complex_seq_analysis/0' num2str(subj, '%02d') '/MVPA/' subfolder '/TDT_results/TDT_results_avg_allreps/' ... 
        'res_accuracy_minus_chance.nii.gz'], 'mask',[base_folder '/Complex_seq_analysis/0' num2str(subj, '%02d') '/MVPA/' subfolder '/mask.nii']);

perm_data = {}; % should be 1 x nperms cell matrix
for i = 1:nperms
    perm_data{i} = cosmo_fmri_dataset([base_folder '/Complex_seq_analysis/0' num2str(subj, '%02d') '/MVPA/' subfolder '/perm/perm' num2str(i, '%03d') '/results_avg_allreps/res_accuracy_minus_chance.nii.gz'], ...
        'mask', [base_folder '/Complex_seq_analysis/0' num2str(subj, '%02d') '/MVPA/' subfolder '/mask.nii']);
end

nbrhood = cosmo_cluster_neighborhood(true_data_nii);
sig_clusters_ds = cosmo_montecarlo_cluster_stat(true_data_nii,nbrhood,'h0_mean',0,'feature_stat', 'none','cluster_stat','tfce','dh',dh,'null',perm_data,'nproc',6);

cosmo_map2fmri(sig_clusters_ds, [base_folder '/Complex_seq_analysis/0' num2str(subj, '%02d') '/MVPA/' subfolder '/sig_clusters_results.nii.gz']);


end
