clear variables;
close all;

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    datalocation = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/trust data recent/scan_behavior/');
    else
        datalocation = glob('?');
    end
end

%% chose models to fit
modelnames = {'ushifted_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;

% get ID list
cd(datalocation{1});
files = dir('trust*.mat');
num_of_subjects = length(files);
%ids = cellstr(NaN(num_of_subjects,1));

%% main loop
%posterior = struct([length(datalocation),length(modelnames)]);
%out = struct([]);
L = [];
parfor ct = 1:num_of_subjects
    filename=files(ct).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(strcmp(id, '212857') || strcmp(id,'881106') ||strcmp(id, '881209'))
        %for m=1:length(modelnames)
            %[posterior(ct,m),out(ct,m)] = clock_sceptic_vba(id(ct),modelnames(m),nbasis, multinomial, multisession, fixed_params_across_runs, fit_propspread);
%            [posterior(ct,5),out(ct,5)] = trust_Qlearning(id, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n);
%            [posterior(ct,1),out(ct,1)] = trust_Qlearning(id, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n);
             [posterior, out] = trust_Qlearning_ushifted(id, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
             L(ct) = out.F;
        %end
    end
end
L_name = sprintf('ushifted_L_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_beta0',multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
save(char(L_name), 'L');
