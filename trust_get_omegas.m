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
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
    else
        datalocation = glob('?');
    end
end


%% choose model's parameters
counter = 1;
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;
utility = 1;

% get ID list
cd(datalocation{1});
files = dir(strcat('*',sprintf('counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d_utility%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, utility),'.mat'));

L = [];
omg = [];
for ct = 1:length(files)
    filename=files(ct).name;
    load(filename)
    L(ct) = out.F;
    omg(ct) = posterior.muTheta(2);
end