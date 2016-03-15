clear variables;
close all;

%needs glob.m (in Project Trust Game/scripts/temporal_instrumental_agent folder

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    pe_location = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        pe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
    else
        pe_location = glob('?');
    end
end

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/';

cd(data_dir_str)
files = dir('*feedback_Times.dat');

num_of_subjects = length(files);

%% choose model's parameters
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;

M_name = strcat(pe_location, sprintf('shiftedPEs_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices));
load(char(M_name));

    trial_2=2;
    total_trials = 192;

for index = 1:num_of_subjects    
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    regressors = load(filename);
    pes = M(M(:,1)==str2double(id),:)';
    pes(1,1) = 0; 
    pes_signed = [regressors(:,1:2), pes];
    pes_signed(trial_2:total_trials,3) = pes_signed(trial_2:total_trials,3) - mean(pes_signed(trial_2:total_trials,3));  %mean-centering  
    pes_unsigned = [regressors(:,1:2), abs(pes)];
    pes_unsigned(trial_2:total_trials,3) = pes_unsigned(trial_2:total_trials,3) - mean(pes_unsigned(trial_2:total_trials,3));  %mean-centering  
    pes_pos = zeros(size(pes));
    pes_neg = zeros(size(pes));
    pes_pos(pes(:,1)>=0,:) = pes(pes(:,1)>=0,:);
    pes_pos(trial_2:total_trials) = pes_pos(trial_2:total_trials) - mean(pes_pos(trial_2:total_trials)); %mean-centering
    pes_neg(pes(:,1)<=0,:) = pes(pes(:,1)<=0,:);
    pes_neg(trial_2:total_trials) = pes_neg(trial_2:total_trials) - mean(pes_neg(trial_2:total_trials)); %mean-centering
    dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'signedPEs' '.dat'],pes_signed,'delimiter','\t','precision','%.6f');
    dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'unsignedPEs' '.dat'],pes_unsigned,'delimiter','\t','precision','%.6f');
    dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'posPEs' '.dat'],[regressors(:,1:2), pes_pos],'delimiter','\t','precision','%.6f');   
    dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'negPEs' '.dat'],[regressors(:,1:2), pes_neg],'delimiter','\t','precision','%.6f');
    dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'pos_neg_PEs' '.dat'],[regressors(:,1:2), pes_pos, pes_neg],'delimiter','\t','precision','%.6f');   
end

% %collinearity
% files = dir('*pos_neg_PEs.dat');
% num_of_subjects = length(files);
% 
% for index = 1:num_of_subjects
%     filename = files(index).name;
%     fprintf('File processing: %s\n', filename);
%     id = filename(isstrprop(filename,'digit'));
%     regressors = load(filename);
%     [R, P] = corrcoef(regressors(:,3:4));
%     corr_matrix(index) = R(1,2);
% end


