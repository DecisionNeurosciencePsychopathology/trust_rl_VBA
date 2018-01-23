

root_dir = 'E:\trust_model_comparision\trust_rl_VBA\';
%dirs = {'all_behav', 'all_scans', 'hc_scan', 'hc_behav'};
dirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};
%subdir = {'all_ages'};
subdir = {};
%specifiy_models = [2 5];
l_struct = [];
for i = 1:length(dirs)
    idx=[];
    files = [];
    files = glob([root_dir dirs{i} filesep subdir{:} filesep '*.mat']);
    for j = 1:length(files)
    %for j = [2 5]
       tmp = load(files{j}); 
       %l_struct.(dirs{i}).L(j,:) = tmp.L;
       l_struct.(dirs{i}).L(j,:) = tmp.out.F;
       idx = [idx j];
    end
    [l_struct.(dirs{i}).posterior,l_struct.(dirs{i}).out] = VBA_groupBMC(l_struct.(dirs{i}).L(idx,:));
    sprintf('%s bor is %d',dirs{i},l_struct.(dirs{i}).out.bor)
    %sprintf('%s ep is %d',dirs{i},l_struct.(dirs{i}).out.ep)
end


