function [] = subject_stats_roi(subfolder, subj_outdir, opt)

arguments
    subfolder
    subj_outdir
    opt.nNulls_per_subject = 500
    opt.subject_num = 'all'
end

% subfolder: folder under Complex_seq_analysis/XXX/MVPA to look for perm and results folders
% subj_outdir: Directory under "subfolder" to save results in (i.e. 'BA44+45')
% nNulls_per_subject: how many nulls from each subject are available to randomly choose from to build group map (default: 500)
% subject_num: subject number or vector of subject numbers (e.g. [1, 2, 3], or 20, default: 'all')

%% sample call:
%% subject_stats_roi('breakpoint_rts_unsmoothed_native_antsmask', 'BA44', 'subject_num', [1, 2, 3])

if exist('/N/slate/brainevo/Implicit_Learning/CoSMoMVPA/mvpa') == 7 % slate
    base_folder = '/N/slate/brainevo/Implicit_Learning';
    addpath(genpath('/N/slate/brainevo/Implicit_Learning/CoSMoMVPA/mvpa'));
    addpath(genpath('/N/slate/brainevo/Implicit_Learning/matlab_nifti_tools'));
else
    base_folder = '/Volumes/data1';
    addpath(genpath('/Users/brainevo/Downloads/CoSMoMVPA/mvpa')); % habilis
    addpath(genpath('/Users/brainevo/Downloads/matlab_nifti_tools'));
end

%%%%%%
%base_folder = '/N/project/brainevo'; % if on slate projects
%%%%%%

% ---------------------------
% Parameters (edit as needed)
% ---------------------------
subj_list = [1 2 3 4 5 6 7 8 9 10 13 14 15 16 17 18 19 20 21 22 23 24];
% Could modify this
% subj_list = [2 4 14 17 18 19]; % good subs
nSubjects = numel(subj_list);
nNulls_per_subject = opt.nNulls_per_subject;
subject_num = opt.subject_num;

if strcmp(subject_num, 'all') == 1
    subject_num = subj_list; % all subs in subj_list
end

for sub = subject_num
    sub = sprintf('%02d', sub);
    res = load([base_folder '/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/' subj_outdir '/TDT_results/TDT_results_avg_allreps/res_accuracy_minus_chance.mat']);
    subj_accuracy = res.avg_all_reps;
    null_vals = zeros(nNulls_per_subject,1);
    for i = 1:nNulls_per_subject
        res = load([base_folder '/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/' subj_outdir '/perm/perm' num2str(i,'%.3d') '/results_avg_allreps/res_accuracy_minus_chance.mat']);
        null_vals(i) = res.avg_all_reps;
    end
    p_val = (sum(null_vals >= subj_accuracy) + 1) / (nNulls_per_subject + 1);
    save([base_folder '/Complex_seq_analysis/0' sub '/MVPA/' subfolder '/' subj_outdir '/subject_accuracy_p-val.mat'], 'p_val');
    fprintf('Done. Subject %s, ROI %s decodes above null at p<%f significance level.\n', sub, subj_outdir, p_val);
end


end
