function compare_all_new_models

root_dir = 'E:\trust_model_comparision\trust_rl_VBA\';
%dirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_corrected', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};
dirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_trustee'};
%dirs={'f_trust_Qlearn_counter_hybrid'};
subdir = {};
%specifiy_models = [2 5];
l_struct = [];

specific_ids = load('trust_hc_subjs_age_filtered');

for i = 1:length(dirs)
    idx=[];
    files = [];
    files = glob([root_dir dirs{i} filesep subdir{:} filesep '*.mat']);
    
    if ~isempty(specific_ids)
        if ~isnumeric(specific_ids)
            %Until I can figure out a better wway for loading in speicific vars
            var_name=fieldnames(specific_ids);
            specific_ids = specific_ids.(var_name{:});
        end
       
       %Create expression string from ids
       specific_ids_str = create_regexp_string(specific_ids);
       
       %Filter files
       file_indices=regexp(files,specific_ids_str);
       
       %Get indices
       index = false(1, numel(file_indices));
       for k = 1:numel(file_indices)
           index(k) = ~isempty(file_indices{k});
       end
       
       files = files(index); %Update files list
       
    end
    
    for j = 1:length(files)
    %for j = [2 5]
       tmp = load(files{j}); 
       %l_struct.(dirs{i}).L(j,:) = tmp.L;
       l_struct.(dirs{i}).L(j,:) = tmp.out.F;
       idx = [idx j];
    end
    %sprintf('%s ep is %d',dirs{i},l_struct.(dirs{i}).out.ep)
end
f_names = fieldnames(l_struct);
l_struct.L_all=[];
for f_name = f_names'
    l_struct.L_all=[l_struct.L_all; l_struct.(f_name{:}).L'];
end

[l_struct.posterior,l_struct.out] = VBA_groupBMC(l_struct.L_all);
sprintf('%s bor is %d',l_struct.out.bor)


function out=create_regexp_string(specific_ids)
out='';
for i = 1:length(specific_ids)
    out = strcat(out,[num2str(specific_ids(i)) '|']);
end
out = out(1:end-1); %Remove trailing or sign
    
    
