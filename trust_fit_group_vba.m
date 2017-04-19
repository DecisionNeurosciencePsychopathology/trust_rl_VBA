clear variables;
close all;

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN') %If using pc
    datalocation = 'E:/Box Sync/Project Trust Game/data/processed/scan_behavior/'; %Set path pointing to single subjects' .mat files
else
    [~, me] = system('whoami'); %If using mac or linux
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/'; %Set path pointing to single subjects' .mat files
    else
        datalocation = '?';
    end
end

%% chose models to fit
modelnames = {'ushifted_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
counter = 1;                    %using counterfactual feedback
multisession = 1;               %modelling runs separately
fixed_params_across_runs = 1;
sigma_kappa = 1;                %kappa (or action bias) parameter
reputation_sensitive = 0;       %modelling trustees' reputation
humanity = 0;                   %modelling humanity
valence_p = 0;                  %modelling valence of trustees
valence_n = 0;
assymetry_choices = 0;          %modelling assymetry in choices
regret = 0;

% get ID list
files = dir([datalocation 'trust*.mat']);
num_of_subjects = length(files);

%% main loop
L = [];
parfor ct = 1:num_of_subjects
    filename=files(ct).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(strcmp(id, '881209')|| strcmp(id, '217909'))
        [posterior, out] = trust_Qlearning_ushifted(datalocation, id, counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
        L(ct) = out.F;
    end
end
L_name = sprintf('L_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
save(char(L_name), 'L');
