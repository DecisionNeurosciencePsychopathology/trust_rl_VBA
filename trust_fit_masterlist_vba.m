clear variables;
close all;

scanner_subjs = 1; %Set to 1 is scanning subjects processed, 0 is behavioral
hallquist_override = 0; %If using Mike's dataset
save_str = ''; %default

os = computer;
if strcmp(os(1:end-2),'PCWIN')
    path_to_trust_ids = 'C:/kod/trust_basic_scripts\';
    if scanner_subjs
        datalocation = glob('E:/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        masterlist = load([path_to_trust_ids 'trust_ids']); %Need a way to automatically update these text files!
    else
        datalocation = glob('E:/Box Sync/Project Trust Game/data/processed/beha_behavior/');
        masterlist = load([path_to_trust_ids 'trust_ids_behav']); %Need a way to automatically update these text files!
        save_str = 'final_behav';
    end
else
    [~, me] = system('whoami');
    me = strtrim(me);
    path_to_trust_ids = '/Users/polinavanyukov/Documents/scripts/Trust/trust_rl_VBA/';
    if strcmp(me,'polinavanyukov')==1
%         datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/');
%         masterlist = load([path_to_trust_ids 'trust_ids_behav']);
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        masterlist = load([path_to_trust_ids 'trust_ids']);
%         masterlist = load([path_to_trust_ids 'trust_hc_subjs_age_filtered']);
    else
        datalocation = glob('?');
    end
end
%masterlist = dlmread([datalocation{:} 'allscansubjs_all_ages.txt']); %Need a way to automatically update these text files!

if hallquist_override
    %path to trust ids is the same currently
    datalocation = glob('C:/Users/wilsonj3/Desktop/hallquist_trust/new_data_1_13/');
    masterlist = load([path_to_trust_ids 'hallquist_trust_ids']); %Need a way to automatically update these text files!
    save_str = 'hallquist/';
end


% get ID list
%cd(datalocation{1});
%files = dir('trust*.mat');
%all_ids = arrayfun(@(x) (x.name(isstrprop(x.name,'digit'))), files, 'UniformOutput', false);

%% chose models to fit
modelnames = {'ushifted_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
sigma_kappa = 1;                %kappa (or action bias) parameter
counter = 0;                    %counter = 1: using counterfactual feedback; counter = 2: using a mixed form of counterfactual+actual rewards
multisession = 0;               %modelling runs separately
fixed_params_across_runs = 1;   
reputation_sensitive = 0;       %modelling trustees' reputation
humanity = 0;                   %modelling humanity
valence_p = 0;                  %modelling valence of trustees
valence_n = 0;                  
assymetry_choices = 0;          %modelling assymetry in choices
regret = 0;

for index=1:length(masterlist.trust_ids)
%for index=1:length(masterlist)
    %id = num2str(masterlist(index));
    id = masterlist.trust_ids(index);
    [posterior, out] = trust_Qlearning_ushifted(datalocation{:}, id, counter, multisession,...
        fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret,save_str);  
    L(index) = out.F;
end

L_name = sprintf('L_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
save(char(L_name), 'L');
