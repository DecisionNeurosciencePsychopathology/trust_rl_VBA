clear variables;
close all;

%needs glob.m (in Project Trust Game/scripts/temporal_instrumental_agent folder

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    values_location = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        pe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        values_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/models mat files/');
        regs_location= glob('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/');
        write_location=strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/',date);
    else
        pe_location = glob('?');
        values_location = glob('?');
        regs_location= glob('?');
        write_location=strcat('?',date);
    end
end

if not(exist(write_location, 'dir'))
    mkdir(write_location);
end


%% choose model's parameters
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;
total_trials = 192;

N_name = strcat(values_location, sprintf('values_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices));
load(char(N_name));

cd(regs_location{1})
files = dir('*feedback_Times.dat');

num_of_subjects = length(files);



for index = 1:num_of_subjects    
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==212857||str2double(id)==215211||str2double(id)==217909 ...
           ||str2double(id)==219471||str2double(id)==881106||str2double(id)==881209)
        if str2double(id) == 46069
            subject = load(strcat(pe_location{1},'trust0',id));
        else
            subject = load(strcat(pe_location{1},'trust',id));
        end
        kept = subject.b.decisions==-1;
        regressors = load(filename);
        values = N(N(:,1)==str2double(id),:)';
        if isempty(N(N(:,1)==str2double(id),:))
            dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/error' id '.dat'], values,'delimiter','\t','precision','%.6f');
        end
        values = values(2:total_trials+1);
        values = circshift(values, total_trials-1);
        values(kept) = values(kept)*-1;
        values(end) = -999;
        values = values - mean(values(1:end-1)); 
        values = [regressors(:,1:2), values];
        dlmwrite([write_location '/trust' id 'values.dat'],values,'delimiter','\t','precision','%.6f');
    end
end
