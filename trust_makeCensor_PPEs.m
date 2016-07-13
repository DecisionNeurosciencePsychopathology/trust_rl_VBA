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
        oldpe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data/');
        pe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        matf_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        regs_location= glob('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/');
        write_location=strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/');
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
regret = 0;
total_trials = 192;

new_M_name = strcat(pe_location,sprintf('modelPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',...
    counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat');
newpes = load(char(new_M_name));
grpnewpes = newpes.M;

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('trust*.mat');
num_of_subjects = length(files);

%hard coded times
dur_choice_display = 300;
dur_feedback = 1200; %1100msec in the file?

%scanning parameters
scan_tr = 1.67;
block_length = 155;
block_end = block_length*scan_tr*1000; %in msec
hemoir = spm_hrf(scan_tr, [6,16,1,1,6,0,32]); % better than resampling and smoothing 
frequency_scale_hz = 10;
% this scale is in msec, but it is separated into bins of X
% Hz (defined by 'frequency_scale' above. the resulting
% output will be in the scale of X Hz.
bin_size = 1/frequency_scale_hz*1000; % convert Hz to ms

for index = 11
%for index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017 ...
        ||str2double(id)==210381||str2double(id)==211669||str2double(id)==216806 ...
        ||str2double(id)==881100||str2double(id)==220024 ...
        ||str2double(id)==220244||str2double(id)==881209||str2double(id)==208572)
        subject = load(strcat(matf_location{1},'trust',id));
        b.regs = [];    
       if not(isfield(subject.b, 'decisions'))
            fprintf('File missing decisions: %s\n', id);
            %Decisions share = 1; keep = -1; no response = 0;
            share =~cellfun(@isempty,strfind(subject.b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(subject.b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(subject.b.PartDecides,'noresponse'));
            subject.b.decisions = zeros(192, 1);
            subject.b.decisions(share) = 1;
            subject.b.decisions(keep) = -1;
            subject.b.decisions(noresponse) = 0;
        end
        kept = subject.b.decisions==-1;
        %% Subject's PEs loaded
       pes = grpnewpes(grpnewpes(:,1)==str2double(id),:)'; 
       pes(1) = 0;
       pes = circshift(pes, length(pes)-2);
       pes = pes(1:192);
       pes(kept) = pes(kept)*-1;       
       PPES_rows = pes > 0; 
       notPPES_rows = pes < 0;
       
       %% PEs decision-aligned
       pes_decision = pes;
       pes_decision = circshift(pes,1);
       pes_decision(1) = 0;
       %plot(pes(1:48)); hold on; plot(pes_decision(1:48),'r');
       ppes_rows_dec = pes_decision > 0;
       ppes_decision = pes_decision;
       ppes_decision(not(ppes_rows_dec)) = 0;
       
       %% GB rows
       GB_rows = strcmp(subject.b.identity,'good') | strcmp(subject.b.identity,'bad');
       notGB_rows = ~GB_rows;
              
        if length(subject.b.ITIfixation_OnsetTime) < 192
            subject.b.ITIfixation_OnsetTime(length(subject.b.ITIfixation_OnsetTime)+1:192)=-999;
            subject.b.ITIfixation_OffsetTime(145:192)=num2cell(-999);
        end
        
        if subject.b.ITIfixation_OnsetTime(192) == -999 || not(iscell(subject.b.firstFixation_OnsetTime))
            firstfix_Onset = subject.b.firstFixation_OnsetTime(1);
            trial2_ITI = 1;
            trial48_ITI = 47;
        else firstfix_Onset = subject.b.ITIfixation_OnsetTime(1);
            trial2_ITI = 2;
            trial48_ITI = 48;
        end

        %for those participants for whom the partnerchoice offset time was not
        %recorded by e-prime
        if iscell(subject.b.partnerchoice_OffsetTime)
            partnerchoice_OffsetTime=subject.b.displaychoice_OnsetTime-1;
        else
            partnerchoice_OffsetTime=subject.b.partnerchoice_OffsetTime;
        end

        fixations = [];
        fixations.event_beg=zeros(48,4);
        fixations.event_end=zeros(48,4); 
        taskness.event_beg =zeros(48,4);
        taskness.event_end=zeros(48,4);
        decision.event_beg=zeros(48,4);
        decision.event_end=zeros(48,4);
        feedback.event_beg=zeros(48,4);
        feedback.event_end=zeros(48,4);
        missed_trials = (subject.b.partnerchoice_RESP==-999);
        notmissed_trials = (missed_trials-1)*-1; 
        
        missed_n_notPPEs = missed_trials|notPPES_rows;
        missed_n_notPPEs_notGB = missed_n_notPPEs|notGB_rows;
        not_miss_n_PPEs = notmissed_trials.*PPES_rows;
        not_miss_n_PPEs_GB = notmissed_trials.*PPES_rows.*GB_rows;
        
        missed_OR_notPPEs_dec = missed_trials|not(ppes_rows_dec);
        not_miss_AND_ppes = not(missed_trials|not(ppes_rows_dec));
        
        %plot(b.notmissed_trials);
        trial1_index = 1;
        trial48_index = 48;
        
        for block= 1:4
            %for fixation screens
            if block == 1 && (subject.b.ITIfixation_OnsetTime(192) == -999 || not(isempty([subject.b.firstFixation_OnsetTime])))
                fixations.event_beg(:,block) = [0; subject.b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset];
            else
                %fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset+block_end*(block-1);
                fixations.event_beg(:,block) = subject.b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset;
            end        
            %fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+block_end*(block-1);
            fixations.event_end(:,block) = subject.b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset;

            %for trial onset to offset; Taskness
            taskness.event_beg(:,block) = subject.b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
            taskness.event_end(:,block) = subject.b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset;

            %for decision onset to response (motor response)
            decision.event_beg(:,block) = subject.b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
            decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 

            %for feedback onset to offset
            feedback.event_beg(:,block) = subject.b.outcome_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
            feedback.event_end(:,block) = subject.b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 
            
            decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 

            %epoch window + missed trials + to censor regressor 
            epoch_window = 0:bin_size:taskness.event_end(48, block);
            
            %% to censor regressors based on the sign of PEs
            %% modify as in main script to include ITIs          
%                 tmp_reg.(['regressors' num2str(block)]).to_censor_PPE = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, not_miss_n_PPEs(trial1_index:trial48_index)); 
%                 tmp_reg.(['regressors' num2str(block)]).to_censor_PPE_GB = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, not_miss_n_PPEs_n_GB(trial1_index:trial48_index));  
            %% modified by AD to include ITI           
            %for PPEs not decision aligned 
            %logical_test = sum(missed_n_notPPEs(trial1_index:trial48_index));
            %include_trials = missed_n_notPPEs(trial1_index:trial48_index);
            
            %for PPEs decision aligned
            logical_test = sum(missed_OR_notPPEs_dec(trial1_index:trial48_index));
            include_trials = missed_OR_notPPEs_dec(trial1_index:trial48_index);
            
            if(logical_test) > 0
                % write 1s for missed trials and not PPEs to censor
                tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, include_trials);
                %tmp_reg.(['regressors' num2str(block)]).to_censor_PPE_GB = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, missed_n_notPPEs_notGB(trial1_index:trial48_index));
            else
                % write a vector of 0s the size of regressors
                tmp_reg.(['regressors' num2str(block)]).to_censor = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, not(include_trials))));
                %tmp_reg.(['regressors' num2str(block)]).to_censor_PPE_GB = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, not_miss_n_PPEs_GB(trial1_index:trial48_index))));
            end
            % flip to_censor
            tmp_reg.(['regressors' num2str(block)]).to_censor = 1-tmp_reg.(['regressors' num2str(block)]).to_censor;
            %tmp_reg.(['regressors' num2str(block)]).to_censor_PPE_GB = 1-tmp_reg.(['regressors' num2str(block)]).to_censor_PPE_GB;
            %% Plotting censor
%             figure(1);
%             subplot(6,1,1);
%             trial_durations = taskness.event_end(:,1)-taskness.event_beg(:,1);
%             trial_durations_10hz = zeros(1,ceil(taskness.event_end(end,1)/100));
%             for epoch=1:length(taskness.event_beg(:,1))
%                 e_begin = taskness.event_beg(epoch,1)/100;
%                 e_end = taskness.event_end(epoch,1)/100;
%                 trial_durations_10hz(e_begin:e_end) = 1;
%             end
%             plot(trial_durations_10hz);
%             
%             subplot(6,1,2);
%             trial_durations_scan = gsresample( ...
%              [zeros(50,1)' trial_durations_10hz], ...
%              10,1./scan_tr);
%             trial_durations_scan = ceil(trial_durations_scan);           
%             trial_durations_scan = [trial_durations_scan ones(1, (block_length-1)-length(trial_durations_scan))];
%             trial_durations_scan = [trial_durations_scan zeros(1,155-length(trial_durations_scan))];
%             plot(trial_durations_scan);
%             
%             subplot(6,1,3);
%             plot(PPES_rows(1:48),'r');
%             %% GSR resample
%             subplot(6,1,4);            
%             plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
%               
%             duration = length(tmp_reg.(['regressors' num2str(block)]).to_censor)*10;
%             vector_sample=1:167:duration;
%             
%             subplot(6,1,5);
%             tmp2 = resample(tmp_reg.(['regressors' num2str(block)]).to_censor,...
%                     10,1);
%             resampled_tocensor = round(tmp2(vector_sample));
%             plot(resampled_tocensor);
%             
%             subplot(6,1,6);
            tmp = gsresample( ...
             [zeros(50,1)' tmp_reg.(['regressors' num2str(block)]).to_censor(1:end-51)], ...
             10,1./scan_tr);
            tmp = round(tmp);
            tmp = [tmp ones(1, (block_length-1)-length(tmp))];
            tmp = [tmp zeros(1,155-length(tmp))];
            tmp_reg.(['regressors' num2str(block)]).to_censor = tmp;
            plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %% Censoring blocks w/ movement
            tmp_reg=censorMovement(id, tmp_reg, block); 
            
            %for next loop iteration, reinitilize variables
            if block < 4
                firstfix_Onset = subject.b.ITIfixation_OnsetTime(trial2_ITI-1+48);
            end
            trial2_ITI=trial2_ITI+48;
            trial48_ITI=trial48_ITI+48;
            trial1_index = trial1_index+48;
            trial48_index = trial48_index+48;  
        end
        to_censor_PPE = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor];
        %to_censor_PPE_GB = [tmp_reg.regressors1.to_censor_PPE_GB tmp_reg.regressors2.to_censor_PPE_GB tmp_reg.regressors3.to_censor_PPE_GB tmp_reg.regressors4.to_censor_PPE_GB];
        to_censor_PPE = transpose(to_censor_PPE);
        %to_censor_PPE_GB = transpose(to_censor_PPE_GB);
        gdlmwrite(strcat(write_location,num2str(id),'to_censor_PPEs_dec'),to_censor_PPE,'\t');
        %gdlmwrite(strcat(write_location,num2str(id),'to_censor_PPEs_GB'),to_censor_PPE_GB,'\t');
    end
end