function just_extract_model_pes_condensed(subdir)

close all;
%subdir={'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT'};
subdir=subdir{:};

%Do we want to zscore the pes?
do_zscore=1;

%subdir = 'f_trust_Qlearn_counter_hybrid';
%If using the null model (f_trust_Qlearn1) do not flip the PEs
if strcmp(subdir,'f_trust_Qlearn1') || strcmp(subdir,'f_trust_Qlearn_null_pmv') || strcmp(subdir,'new_f_trust_Qlearn_null_pmv')
    flipthis=0;
else
    flipthis=1;
end

    pe_location = glob(['E:\trust_model_comparision\trust_rl_VBA\' subdir '\PEs']);
    matf_location = glob('E:\Box Sync\Project Trust Game\data\processed\scan_behavior\');
    regs_location= glob('E:\data\trust\regs\');
    write_location=strcat(['E:\trust_model_comparision\trust_rl_VBA\' subdir '\regs\'],date);
if not(exist(write_location, 'dir'))
    mkdir(write_location);
end

if strcmp(subdir,'f_trust_Qlearn1') ||  strcmp(subdir,'f_trust_Qlearn_counter_hybrid_regret_pmv') || strcmp(subdir,'f_trust_Qlearn_null_pmv')%Band-aid, but really it should pull the params from the file name
    counter = 0;
else
    counter = 1;
end

kept_hist = [];
pes_flipped_hist = [];
model_str = glob([pe_location{:} 'model*.mat']);
newpes=load(model_str{:});
%newpes=load('E:\trust_model_comparision\trust_rl_VBA\f_trust_Qlearn_policy_censor0_mltrun1_kS_kT\PEs\modelPEs_f_trust_Qlearn_policy_censor0_mltrun1_kS_kT.mat');
grpnewpes = newpes.M;

cd(regs_location{1})
files = dir('*feedback_Times.dat');
num_of_subjects = length(files);

for index = 1:num_of_subjects
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==202021||str2double(id)==219956||str2double(id)==220017||...
            str2double(id)==210381||str2double(id)==211669||str2double(id)==216806||...
            str2double(id)==881100||str2double(id)==220024||...
            str2double(id)==220244||str2double(id)==881209)
        if str2double(id) == 46069
            %subject = load(strcat(matf_location{1},'trust0',id));
        else
            subject = load(strcat(matf_location{1},'trust',id));
        end
        %check if decisions field has been created as part of the subject's
        %struct data, if not, create one
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
        
        %% trustee contrasts vectors
        trustee_GB = zeros(length(subject.b.identity),1);
        trustee_G0 = zeros(length(subject.b.identity),1);
        trustee_HC = ones(length(subject.b.identity),1);
        %trustee_BA = ones(length(subject.b.identity),1)*-1;
        trustee_BA = zeros(length(subject.b.identity),1);
        trustee_action = ones(length(subject.b.identity),1);
        trustee_survey_GB = zeros(length(subject.b.identity),1);
        
        %% human vs non-human
        trustee_HC(strcmp(subject.b.identity,'computer')) = 0;
        
        %% good and bad
        trustee_GB(strcmp(subject.b.identity,'good')) = 1;
        trustee_GB(strcmp(subject.b.identity,'bad')) = -1;
        trustee_BA(strcmp(subject.b.identity,'bad')) = 1;
        trustee_G0(strcmp(subject.b.identity,'good')) = 1;
        
        %% Continuantion of PE calc
        kept = subject.b.decisions==-1;
        regressors = load(filename);
        pes = grpnewpes(grpnewpes(:,1)==str2double(id),:)';
        if isempty(pes)
            fprintf([id ' is not in the HC list\n'])
            continue
        end
        kept = subject.b.decisions==-1;
        regressors = load(filename);
        pes(1) = 0;
        kept = subject.b.decisions==-1;
        pes = circshift(pes, length(pes)-2);
        
        %Let's see the difference for models
%         figure(index)
%         clf;
%         plot(pes)
%         hold on
%         pes_flipped = pes;
%         pes_flipped(kept) = pes(kept)*-1;
%         plot(pes_flipped)
%         kept_hist = [kept_hist kept];
%         pes_flipped_hist.(['id' id]) = pes(kept)*-1;
        
        
        if flipthis
            pes(kept) = pes(kept)*-1;
        end
        pes = pes(1:192);
        [pes_r, pes_c] = size(pes);
        if pes_r ~= length(pes)
            pes=pes';
            calcpes=calcpes';
            ppes=ppes';
        end
        
        pes = [regressors(:,1:2), pes];
        
        %10/10/2017 Zscore the pes before they go into the regressors to
        %put betas o nthe same scale
        if do_zscore
            pes(:,3) = zscore(pes(:,3));
        else
            pes(:,3) = pes(:,3) - mean(pes(:,3));
        end
        
        
        dlmwrite([write_location '/trust' id 'PEs' '.dat'], pes,'delimiter','\t','precision','%.6f'); %Original we always use!
        
        
        
        
% %         %% Contrasts
% %         pesxBA = [regressors(:, 1:2), pes(:,3).* trustee_BA];
% %         pesxGB = [regressors(:, 1:2), pes(:,3).* trustee_GB];
% %         pesxG0 = [regressors(:, 1:2), pes(:,3).* trustee_G0];
% %         %pesxHC = [regressors(:, 1:2), pes(:,3).* trustee_HC];
% %         
% %         %Mean center
% %         GB_rows = not(trustee_GB == 0);
% %         pesxBA(:,3) = pesxBA(:,3) - mean(pesxBA(1:length(pes),3));
% %         pesxGB(GB_rows,3) = pesxGB(GB_rows,3) - mean(pesxGB(GB_rows,3));
% %         pesxG0(:,3) = pesxG0(:,3) - mean(pesxG0(:,3));
% %         %pesxHC(:,3) = pesxHC(:,3) - mean(pesxHC(:,3));
% %         
% %         
% %         dlmwrite([write_location '/trust' id 'BAxPEs' '.dat'],pesxBA,'delimiter','\t','precision','%.6f');
% %         dlmwrite([write_location '/trust' id 'GBxPEs' '.dat'],pesxGB,'delimiter','\t','precision','%.6f');
% %         dlmwrite([write_location '/trust' id 'G0xPEs' '.dat'],pesxG0,'delimiter','\t','precision','%.6f');
% %         %dlmwrite([write_location '/trust' id 'HCxPEs' '.dat'],pesxHC,'delimiter','\t','precision','%.6f');
    end
end