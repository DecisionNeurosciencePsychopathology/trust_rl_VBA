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
        regs_location= glob('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/');
        write_location=strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/',date);
    else
        pe_location = glob('?');
        regs_location= glob('?');
        write_location=strcat('?',date);
    end
end

if not(exist(write_location, 'dir'))
    mkdir(write_location);
end
        
    

%% choose model's parameters and load model
counter = 1;
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;
total_trials = 192;

M_name = strcat(pe_location, 'models mat files/',sprintf('PEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices));
load(char(M_name));



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
        %check if decisions field has been created as part of the subject's
        %struct data, if not, create one
        if not(isfield(subject.b, 'decisions'))
            fprintf('File missing decisions: %s\n', id);
            %Decisions share = 1; keep = -1; no response = 0;
            share =~cellfun(@isempty,strfind(subject.b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(subject.b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(subject.b.PartDecides,'noresponse'));
            subject.b.decisions = zeros(trials-start, 1);
            subject.b.decisions(share) = 1;
            subject.b.decisions(keep) = -1;
            subject.b.decisions(noresponse) = 0;
        end
        kept = subject.b.decisions==-1;
        regressors = load(filename);
        pes = M(M(:,1)==str2double(id),:)';
        if isempty(M(M(:,1)==str2double(id),:))
            dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/error' id '.dat'], pes,'delimiter','\t','precision','%.6f');
        end
        pes = circshift(pes, total_trials-1); %shifts the array so that the first element of pe = pe(1) and last element of pe = id
        pes(end)=0;
        pes_signed = [regressors(:,1:2), pes];
        
        pes_actions = pes_signed;
        pes_actions(kept,3) = pes_actions(kept,3)*-1; %flipping PEs for share when the subject actually decided to keep
        
        pes_signed(1:total_trials-1,3) = pes_signed(1:total_trials-1,3) - mean(pes_signed(1:total_trials-1,3));  %mean-centering 
        pes_unsigned = [regressors(:,1:2), abs(pes)];
        pes_unsigned(1:total_trials-1,3) = abs(pes_unsigned(1:total_trials-1,3) - mean(pes_unsigned(1:total_trials-1,3)));  %mean-centering  
        pes_pos = zeros(size(pes));
        pes_neg = zeros(size(pes));
        pes_action_pos = zeros(size(pes_actions));
        pes_action_neg = zeros(size(pes_actions));
        pes_pos(pes(:,1)>=0,:) = pes(pes(:,1)>=0,:);
        pes_pos(1:total_trials-1) = pes_pos(1:total_trials-1) - mean(pes_pos(1:total_trials-1)); %mean-centering
        pes_neg(pes(:,1)<=0,:) = pes(pes(:,1)<=0,:);
        pes_neg(1:total_trials-1) = pes_neg(1:total_trials-1) - mean(pes_neg(1:total_trials-1)); %mean-centering
        pes_action_pos(pes_actions(:,3)>0,3) = pes_actions(pes_actions(:,3)>0,3);
        pes_action_pos(1:total_trials-1,3) = pes_action_pos(1:total_trials-1,3) - mean(pes_action_pos(1:total_trials-1,3)); %mean-centering
        pes_neg(pes(:,1)<=0,:) = pes(pes(:,1)<=0,:);
        pes_neg(1:total_trials-1) = pes_neg(1:total_trials-1) - mean(pes_neg(1:total_trials-1)); %mean-centering
        pes_action_neg(pes_actions(:,3)<=0,3) = pes_actions(pes_actions(:,3)<=0,3);
        pes_action_neg(1:total_trials-1,3) = pes_action_neg(1:total_trials-1,3) - mean(pes_action_neg(1:total_trials-1,3)); %mean-centering
        pes_actions(1:total_trials-1,3) = pes_actions(1:total_trials-1,3) - mean(pes_actions(1:total_trials-1,3));  %mean-centering 
        %dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'signedPEs' '.dat'],pes_signed,'delimiter','\t','precision','%.6f');
        dlmwrite([write_location '/trust' id 'unsignedPEs' '.dat'],pes_unsigned,'delimiter','\t','precision','%.6f');
        dlmwrite([write_location '/trust' id 'posPEs' '.dat'],[regressors(:,1:2), pes_pos],'delimiter','\t','precision','%.6f');   
        dlmwrite([write_location '/trust' id 'negPEs' '.dat'],[regressors(:,1:2), pes_neg],'delimiter','\t','precision','%.6f');
        %dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted' id 'pos_neg_PEs' '.dat'],[regressors(:,1:2), pes_pos, pes_neg],'delimiter','\t','precision','%.6f');
        %plot(pes_actions(:,3));
        dlmwrite([write_location '/trust' id 'action_PEs' '.dat'], pes_actions,'delimiter','\t','precision','%.6f');
        %dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/trust_shifted'...
        %    id 'pos_neg_action_PEs' '.dat'],[regressors(:,1:2), pes_action_pos(:,3), pes_action_neg(:,3)],'delimiter','\t','precision','%.6f');
        dlmwrite([write_location '/trust'...
             id 'pos_action_PEs' '.dat'],[regressors(:,1:2), pes_action_pos(:,3)],'delimiter','\t','precision','%.6f');
        dlmwrite([write_location '/trust'...
        id 'negs_action_PEs' '.dat'],[regressors(:,1:2), pes_action_neg(:,3)],'delimiter','\t','precision','%.6f');
    end
end

% %collinearity
% files = dir('*pos_neg_action_PEs.dat');
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
% 

