function trust_extract_PEs(subdirs)
close all;

%subdirs = {'f_trust_Qlearn_counter_hybrid' 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_trustee'};
%subdirs = {'f_trust_Qlearn_counter_hybrid'};
%subdirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_corrected', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};

%Really need to fix this cd stuff
current_dir=pwd;


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
    if strcmp(subdirs{j},'f_trust_Qlearn1') ||  strcmp(subdirs{j},'f_trust_Qlearn_counter_hybrid_regret_pmv') || strcmp(subdirs{j},'f_trust_Qlearn_null_pmv')%Band-aid, but really it should pull the params from the file name
        counter = 0;
    else
        counter = 1;
    end
    multisession = 1;
    fixed_params_across_runs = 1;
    sigma_kappa = 1;
    reputation_sensitive = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
    regret = 0;
    
    file_str = sprintf('cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    
    % get ID list
    cd(datalocation{1});
    %files = dir(strcat('*',sprintf('counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d_regret%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat'));
    files = dir(strcat('*',file_str,'.mat'));
    
    %Rename files to where they came from if long form fails
    if isempty(files)
        files = dir(strcat('*.mat'));
        file_str = subdirs{j};
    end
    
    num_of_subjects = length(files);
%     P=zeros(num_of_subjects,193); %CALCULATED PEs
%     M=zeros(num_of_subjects,193); %MODEL PEs
%     N=zeros(num_of_subjects,193); %value
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
        
        if multisession
            try
                nt = length(out.suffStat.muX(2,:))/blocklength;
                count=2;
                M_tmp=[];
                for j=1:nt
                   M_tmp = [M_tmp out.suffStat.muX(count,blocklength*(j-1)+1:(blocklength*j))];
                   count = count+2;
                end
                 M(ct,:) = [subject_id,M_tmp];
            catch
                warning('Probably missing trials or blocks for subject %d',subject_id)
                zeros_to_add=length(M)-length(M_tmp);
                M_tmp = [M_tmp zeros(1,zeros_to_add-1)];
                M(ct,:) = [subject_id,M_tmp];
            end
            
        end
        
        
    end
    
    if ~exist('PEs', 'dir')
        mkdir('PEs')
    end
    
    
    %P_name = ['PEs' filesep sprintf('calcPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    P_name = ['PEs' filesep 'calcPEs_' file_str];
    save(char(P_name), 'P');
    
   % M_name = sprintf('modelPEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    %M_name = ['PEs' filesep sprintf('modelPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    M_name = ['PEs' filesep 'modelPEs_' file_str];
    save(char(M_name), 'M');
    
    %N_name = sprintf('values_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
    %N_name = ['PEs' filesep sprintf('values_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret)];
    N_name = ['PEs' filesep 'values_' file_str];
    save(char(N_name), 'N');
    close all;
end

cd(current_dir)