clear variables;
close all;

os = computer;
if strcmp(os(1:end-2),'PCWIN')
    datalocation = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
    else
        datalocation = glob('?');
    end
end

masterlist = dlmread('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/new subjects.txt');

% get ID list
%cd(datalocation{1});
%files = dir('trust*.mat');
%all_ids = arrayfun(@(x) (x.name(isstrprop(x.name,'digit'))), files, 'UniformOutput', false);

%% chose models to fit
modelnames = {'ushifted_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
counter = 1;                    %using counterfactual feedback
multisession = 0;               %modelling runs separately
fixed_params_across_runs = 1;   
sigma_kappa = 1;                %kappa (or action bias) parameter
reputation_sensitive = 0;       %modelling trustees' reputation
humanity = 0;                   %modelling humanity
valence_p = 0;                  %modelling valence of trustees
valence_n = 0;                  
assymetry_choices = 0;          %modelling assymetry in choices
regret = 0;


parfor index=1:length(masterlist)
%for index=1:length(masterlist)
    id = num2str(masterlist(index));
    [posterior, out] = trust_Qlearning_ushifted(id, counter, multisession,...
        fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);  
    L(index) = out.F;
end

% L_name = sprintf('L_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
% save(char(L_name), 'L');