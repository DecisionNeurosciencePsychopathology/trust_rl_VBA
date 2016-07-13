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
        matf_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
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

cd(matf_location{1})
files = dir('trust*.mat');
num_of_subjects = length(files);

scan_tr = 1.67;
frequency_scale_hz = 10;
block_length = 155;
block_end = block_length*scan_tr*1000; 
% this scale is in msec, but it is separated into bins of X
% Hz (defined by 'frequency_scale' above. the resulting
% output will be in the scale of X Hz.
bin_size = 1/frequency_scale_hz*1000;

%for ind = 15
for ind=1:num_of_subjects
    filename = files(ind).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017 ...
        ||str2double(id)==210381||str2double(id)==211669||str2double(id)==216806 ...
        ||str2double(id)==881100||str2double(id)==220024 ...
        ||str2double(id)==220244||str2double(id)==881209)
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
        
        taskness.event_beg =zeros(48,4);
        taskness.event_end=zeros(48,4);
        if iscell(subject.b.partnerchoice_RESP)
            missed_trials = subject.b.decisions==0;
        else
            missed_trials = (subject.b.partnerchoice_RESP==-999);
        end
        missed_blockbeg = missed_trials;
        missed_blockbeg([1 49 97 145]) = 1;        
        
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
        trial1_index = 1;
        trial48_index = 48;
        
        %% block loop
        for block =1:4            
            %for trial onset to offset; Taskness
            taskness.event_beg(:,block) = subject.b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
            taskness.event_end(:,block) = subject.b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset;
            epoch_window = 0:bin_size:taskness.event_end(48, block);
            
            %% CENSOR vector modified by AD to include ITI
            if(sum(missed_blockbeg(trial1_index:trial48_index))) > 0
                % write 1s for missed trials and not PPEs to censor
                tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, missed_blockbeg(trial1_index:trial48_index));
            else
                % write a vector of 0s the size of regressors
                tmp_reg.(['regressors' num2str(block)]).to_censor = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, not(missed_blockbeg(trial1_index:trial48_index)))));
            end
            %% flip to_censor
            tmp_reg.(['regressors' num2str(block)]).to_censor = 1-tmp_reg.(['regressors' num2str(block)]).to_censor;
            
            %% resampling
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
            
            %% for next loop iteration, reinitilize variables
            if block < 4
                firstfix_Onset = subject.b.ITIfixation_OnsetTime(trial2_ITI-1+48);
            end
            trial2_ITI=trial2_ITI+48;
            trial48_ITI=trial48_ITI+48;
            trial1_index = trial1_index+48;
            trial48_index = trial48_index+48;       
        end
        to_censor_PE_d = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor];
        %to_censor_PPE_GB = [tmp_reg.regressors1.to_censor_PPE_GB tmp_reg.regressors2.to_censor_PPE_GB tmp_reg.regressors3.to_censor_PPE_GB tmp_reg.regressors4.to_censor_PPE_GB];
        to_censor_PE_d = transpose(to_censor_PE_d);
        %to_censor_PPE_GB = transpose(to_censor_PPE_GB);
        gdlmwrite(strcat(write_location,num2str(id),'to_censor_PEs_d'),to_censor_PE_d,'\t');
    end
end