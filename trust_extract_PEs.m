clear variables;
close all;

%subdirs = {'f_trust_Qlearn_counter_hybrid' 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_trustee'};
subdirs = {'f_trust_Qlearn_counter_hybrid'};

for j = 1:length(subdirs)
    %Quick username check, and path setting, this may have to change depending
    %on the machine you are currently working on!
    os = computer;
    if strcmp(os(1:end-2),'PCWIN')
        datalocation = glob(['E:\trust_model_comparision\trust_rl_VBA\' subdirs{j} filesep]);
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
    regret = 0;
    
    % get ID list
    cd(datalocation{1});
    %files = dir(strcat('*',sprintf('counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d_regret%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat'));
    files = dir(strcat('*',sprintf('cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat'));
    num_of_subjects = length(files);
    P=zeros(num_of_subjects,193); %CALCULATED PEs
    M=zeros(num_of_subjects,193); %MODEL PEs
    N=zeros(num_of_subjects,193); %value
    for ct = 1:num_of_subjects
        filename=files(ct).name;
        fprintf('File processing: %s\n', filename);
        %subject_id = filename(isstrprop(filename,'digit'));
        load(filename);
        if ischar(id)
            subject_id = str2double(id);
        else
            subject_id = id;
        end
        P(ct,:) = [subject_id, out.suffStat.PE];
        M(ct,:) = [subject_id, out.suffStat.muX(2,:)];
        N(ct,:) = [subject_id, out.suffStat.muX(1,:)];
        %close all;
    end
    
    if ~exist('PEs', 'dir')
        mkdir('PEs')
    end
    
    
    P_name = ['PEs' filesep sprintf('calcPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    save(char(P_name), 'P');
    
   % M_name = sprintf('modelPEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    M_name = ['PEs' filesep sprintf('modelPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    save(char(M_name), 'M');
    
    %N_name = sprintf('values_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
    N_name = ['PEs' filesep sprintf('values_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    save(char(N_name), 'N');
    close all;
end