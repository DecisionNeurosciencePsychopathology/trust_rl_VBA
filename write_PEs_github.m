clear variables;
close all;


subdir = 'f_trust_Qlearn_counter_hybrid';
flipthis=0;

%needs glob.m (in Project Trust Game/scripts/temporal_instrumental_agent folder

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    %oldpe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data/');
    pe_location = glob(['E:\trust_model_comparision\trust_rl_VBA\' subdir '\PEs']);
    matf_location = glob('E:\Box Sync\Project Trust Game\data\processed\scan_behavior\');
    regs_location= glob('E:\data\trust\regs\feedback_times\');
    write_location=strcat(['E:\trust_model_comparision\trust_rl_VBA\' subdir '\regs\'],date);
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        oldpe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data/');
        pe_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        matf_location = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        regs_location= glob('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/July 17/');
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
regret = 0;
total_trials = 192;
%M_name = strcat(oldpe_location,sprintf('PEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat');
K_name = strcat(pe_location,sprintf('calcPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat');
%oldpes = load(char(M_name));
calcpes = load(char(K_name));
% N_name = strcat(pe_location, sprintf('modelvalues_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',...
%     counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices),'.mat');
N_name = strcat(pe_location, sprintf('values_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',...
    counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices,regret),'.mat');
values = load(char(N_name));

% new_M_name = strcat(pe_location,sprintf('modelPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',...
%     counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat');

new_M_name =  strcat(pe_location,sprintf('modelPEs_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',...
      counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret),'.mat');

newpes = load(char(new_M_name));
grpnewpes = newpes.M;
%grpoldpes = oldpes.M;
grpvalues = values.N;
grpcalcpes = calcpes.P;

cd(regs_location{1})
files = dir('*feedback_Times.dat');
%files = dir('*decision_Times.dat');
num_of_subjects = length(files);

%for index = 60
for index = 1:num_of_subjects
    try
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==202021||str2double(id)==219956||str2double(id)==220017||...
            str2double(id)==210381||str2double(id)==211669||str2double(id)==216806||...
            str2double(id)==881100||str2double(id)==220024||...
            str2double(id)==220244||str2double(id)==881209)
        if str2double(id) == 46069
            subject = load(strcat(matf_location{1},'trust0',id));
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
        kept = subject.b.decisions==-1;
        regressors = load(filename);
        
        %% the PEs are not found for the following subjects
        %if isempty(M(M(:,1)==str2double(id),:))
        if isempty(grpnewpes(grpnewpes(:,1)==str2double(id),:))
            dlmwrite(['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/error' id '.dat'], pes,'delimiter','\t','precision','%.6f');
        end
        %% Subject's PEs loaded
        pes = grpnewpes(grpnewpes(:,1)==str2double(id),:)';
        %oldpes = grpoldpes(grpoldpes(:,1)==str2double(id),:)';
        calcpes= grpcalcpes(grpcalcpes(:,1)==str2double(id),:)';
        values = grpvalues(grpvalues(:,1)==str2double(id),:)';
        
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
        
        %% trustee actions
        trustee_action(strcmp(subject.b.TrusteeDecides,'keep'))=-1;
        
        %% subject's decisions
        % share - keep;
        [rows, columns] = size(subject.b.decisions);
        if columns == 192
            subject.b.decisions = subject.b.decisions';
        end
        shareVSkeep = zeros(size(subject.b.partnerchoice_RESP));
        shareVSkeep(subject.b.decisions==1 & subject.b.partnerchoice_RESP ~= -999) = 1;
        %shareVSkeep(subject.b.decisions==-1 & subject.b.partnerchoice_RESP ~= -999) = -1;
        
        %% survey based vectors (not using)
        %         load('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/surveys/trust_group_survey.mat');
        %         if ischar(id)
        %             rows = group_survey_table.ID==str2double(id);
        %         else
        %             rows = group_survey_table.ID==id;
        %         end
        %         survey = group_survey_table(rows,:);
        %         good_rating = survey.Rating_Pre(strcmp(survey.Trustee,'good') & strcmp(survey.Question,'trustworthy'));
        %         bad_rating = survey.Rating_Pre(strcmp(survey.Trustee,'bad') & strcmp(survey.Question,'trustworthy'));
        %         trustee_survey_GB(strcmp(subject.b.identity,'good')) = 1*good_rating/7;
        %         trustee_survey_GB(strcmp(subject.b.identity,'bad')) = -1*bad_rating/7;
        %% aligning PES
        %getting rid of the ID
        pes(1) = 0;
        %oldpes(1) = 0;
        calcpes(1) =0;
        values(1) = 0;
        %
        %        figure(97); clf;
        %        subplot(3,1,1);
        %        hold on;
        %        %plot(oldpes(1:48),'b');                             %OLD PES
        %        %plot(circshift(oldpes(1:48),total_trials+1),'b');    %OLD PES ALIGNED w/ NEW PES
        %        plot(circshift(oldpes(1:48),total_trials-1),'b');
        %        %plot(calcpes(1:48),'r');                             %CALCULATED PES
        %        plot(pes(1:48),'k');                                 %MODEL-TRACKED PES
        %        subplot(3,1,2);
        %        %plot(values(1:48),'r');                             %MODEL-tracked VALUE
        %        plot(circshift(values(1:48), total_trials+1),'r');   %VALUE ALIGNED W/ DECISIONS MADE
        %        subplot(3,1,3);
        %        hold on;
        %        plot(subject.b.decisions(1:48), 'b');
        %        plot(trustee_action(1:48),'k');
        %        hold off;
        
        %         oldpes_align2new = circshift(oldpes,length(oldpes)+1);   %ALIGNING OLD PES W/ NEW PES
        %         oldpes_align2new(1) = 0;
        
        % aligning new PEs w/ old PEs
        %       oldpes = circshift(oldpes, length(oldpes)-1);
        pes = circshift(pes, length(pes)-2);
        calcpes = circshift(calcpes, length(calcpes)-2);
        values = circshift(values, length(values)-2); %expected value is aligned with the beginning of next trial's feedback
        
        %        corr_oldpes_calcpes = corr(oldpes_align2new, calcpes(1:192));
        %        corr_oldpes_modelpes = corr(oldpes_align2new, pes(1:192));
        
        %       old_old = load(strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/31-Mar-2016/trust', id, 'action_PEs.dat'));
        
        %        figure(96); clf;
        %        subplot(3,1,1);
        %        hold on;
        %        %plot(oldpes(1:48),'b');                             %OLD PES
        %        plot(pes(1:48),'k');                                 %MODEL-TRACKED PES
        %        subplot(3,1,2);
        %        plot(values(1:48),'r');                             %MODEL-tracked VALUE
        %        %plot(circshift(values(1:48), total_trials+1),'r');   %VALUE ALIGNED W/ DECISIONS MADE
        %        subplot(3,1,3);
        %        hold on;
        %        plot(subject.b.decisions(1:48), 'b');
        %        plot(trustee_action(1:48),'k');
        %        hold off;
        %
        %% Flipping all PEs for share when subject actually decided to keep
        if flipthis
            pes(kept) = pes(kept)*-1;
            calcpes(kept) = calcpes(kept)*-1;
        end
        %       oldpes(kept) = oldpes(kept)*-1;
        
        %% value and PE assembly w/ time epochs
        % Shortening PEs to 192 trials
        pes = pes(1:192);
        calcpes = calcpes(1:192);
        values = values(1:192);
        values_abs=abs(values);
        
        %If the pes somwhow got flipped dimmensionally
        [pes_r, pes_c] = size(pes);        
        if pes_r ~= length(pes)
            pes=pes';
            calcpes=calcpes';
            ppes=ppes';
        end
        
        
        %% PE, values and PPEs aligned w/ decision of trial t + 1
        pes_decision = pes;
        pes_decision = circshift(pes,1);
        pes_decision(1) = 0;
        %plot(pes(1:48)); hold on; plot(pes_decision(1:48),'r');
        ppes_rows_dec = pes_decision > 0;
        ppes_decision = pes_decision;
        ppes_decision(not(ppes_rows_dec)) = 0;
        values_dec = values;
        values_dec = circshift(values,1);
        
        %% VALIDATED UP TO THIS POINT WITH SIGNED PE MAPS
        PPES_rows = pes > 0;
        ppes = pes;
        ppes(not(PPES_rows)) = 0;
                
        pes = [regressors(:,1:2), pes];
        calcpes = [regressors(:,1:2), calcpes];
        %oldpes = [regressors(:,1:2), oldpes];
        ppes = [regressors(:,1:2), ppes];
        %ppes2 = [regressors(:,1:2), ppes];
        values = [regressors(:,1:2), values];
        values_abs = [regressors(:,1:2), values_abs];
        values_dec_abs = abs(values_dec);
        
        %% Interaction regressors
        %% 2-way interactions
        pesxBA = [regressors(:, 1:2), pes(:,3).* trustee_BA];
        pesxGB = [regressors(:, 1:2), pes(:,3).* trustee_GB];
        pesxG0 = [regressors(:, 1:2), pes(:,3).* trustee_G0];
        pesxHC = [regressors(:, 1:2), pes(:,3).* trustee_HC];
        pesxTrial = pes(:,3).*(1:192)';
        
        ppesxBA = [regressors(:, 1:2), ppes(:,3).* trustee_BA];
        ppesxGB = [regressors(:, 1:2), ppes(:,3).* trustee_GB];
        ppesxG0 = [regressors(:, 1:2), ppes(:,3).* trustee_G0];
        ppesxHC = [regressors(:, 1:2), ppes(:,3).* trustee_HC];
        
        %         pesxsurveyGB = [regressors(:, 1:2), pes_actions(:,3).* trustee_survey_GB];
        pesxp_decision = [regressors(:, 1:2), pes(:,3).* shareVSkeep];
        ppesxp_decision =  [regressors(:, 1:2), ppes(:,3).* shareVSkeep];
        
        pes_decxBA = pes_decision.* trustee_BA;
        pes_decxGB = pes_decision.* trustee_GB;
        pes_decxG0 = pes_decision.* trustee_G0;
        pes_decxHC = pes_decision.* trustee_HC;
        pes_decxpdecision = pes_decision.*shareVSkeep;
        pes_decxblockOrder = pes_decision.*[ones(48,1); 2*ones(48,1); 3*ones(48,1); 4*ones(48,1)];
        
        ppes_decxBA = ppes_decision.* trustee_BA;
        ppes_decxGB = ppes_decision.* trustee_GB;
        ppes_decxG0 = ppes_decision.* trustee_G0;
        ppes_decxHC = ppes_decision.* trustee_HC;
        ppes_decxpdecision = ppes_decision.*shareVSkeep;
        
        
        %% 3-way interactions
        pesxBAxp_decision =[regressors(:, 1:2), pesxBA(:,3).*shareVSkeep];
        pesxGBxp_decision =[regressors(:, 1:2), pesxGB(:,3).*shareVSkeep];
        pesxG0xp_decision =[regressors(:, 1:2), pesxG0(:,3).*shareVSkeep];
        pesxHCxp_decision =[regressors(:, 1:2), pesxHC(:,3).*shareVSkeep];
        
        ppesxBAxp_decision =[regressors(:, 1:2), ppesxBA(:,3).*shareVSkeep];
        ppesxGBxp_decision =[regressors(:, 1:2), ppesxGB(:,3).*shareVSkeep];
        ppesxG0xp_decision =[regressors(:, 1:2), ppesxG0(:,3).*shareVSkeep];
        ppesxHCxp_decision =[regressors(:, 1:2), ppesxHC(:,3).*shareVSkeep];
        
        pes_decxBAxpdecision = pes_decxBA.*shareVSkeep;
        pes_decxGBxpdecision = pes_decxGB.*shareVSkeep;
        pes_decxG0xpdecision = pes_decxG0.*shareVSkeep;
        pes_decxHCxpdecision = pes_decxHC.*shareVSkeep;
        pes_decxBAxOrder = pes_decxBA.*[ones(48,1); 2*ones(48,1); 3*ones(48,1); 4*ones(48,1)];
        
        ppes_decxBAxpdecision = ppes_decxBA.*shareVSkeep;
        ppes_decxGBxpdecision = ppes_decxGB.*shareVSkeep;
        ppes_decxG0xpdecision = ppes_decxG0.*shareVSkeep;
        ppes_decxHCxpdecision = ppes_decxHC.*shareVSkeep;
        
        %% MEAN_CENTERING
        GB_rows = not(trustee_GB == 0);
        
        pesxBA(:,3) = pesxBA(:,3) - mean(pesxBA(1:length(pes),3));
        pesxGB(GB_rows,3) = pesxGB(GB_rows,3) - mean(pesxGB(GB_rows,3));
        pesxG0(:,3) = pesxG0(:,3) - mean(pesxG0(:,3));
        pesxHC(:,3) = pesxHC(:,3) - mean(pesxHC(:,3));
        pesxTrial = pesxTrial - mean(pesxTrial);
        
        pesxp_decision(:,3) = pesxp_decision(:,3) - mean(pesxp_decision(:,3));
        pesxBAxp_decision(:,3) = pesxBAxp_decision(:,3) - mean(pesxBAxp_decision(:,3));
        pesxGBxp_decision(GB_rows,3) = pesxGBxp_decision(GB_rows,3) - mean(pesxGBxp_decision(GB_rows,3));
        pesxG0xp_decision(:,3) = pesxG0xp_decision(:,3) - mean(pesxG0xp_decision(:,3));
        pesxHCxp_decision(:,3) = pesxHCxp_decision(:,3) - mean(pesxHCxp_decision(:,3));
        
        ppesxBA(PPES_rows,3) = ppesxBA(PPES_rows,3) - mean(ppesxBA(PPES_rows,3));
        
        ppesxGB(PPES_rows,3) = ppesxGB(PPES_rows,3) - mean(ppesxGB(PPES_rows & GB_rows,3));
        ppesxG0(PPES_rows,3) = ppesxG0(PPES_rows,3) - mean(ppesxG0(PPES_rows,3));
        ppesxHC(PPES_rows,3) = ppesxHC(PPES_rows,3) - mean(ppesxHC(PPES_rows,3));
        
        ppesxp_decision(PPES_rows,3) = ppesxp_decision(PPES_rows,3) - mean(ppesxp_decision(PPES_rows,3));
        ppesxBAxp_decision(PPES_rows,3) = ppesxBAxp_decision(PPES_rows,3) - mean(ppesxBAxp_decision(PPES_rows,3));
        ppesxGBxp_decision(PPES_rows,3) = ppesxGBxp_decision(PPES_rows,3) - mean(ppesxGBxp_decision(PPES_rows & GB_rows,3));
        ppesxG0xp_decision(PPES_rows,3) = ppesxG0xp_decision(PPES_rows,3) - mean(ppesxG0xp_decision(PPES_rows,3));
        ppesxHCxp_decision(PPES_rows,3) = ppesxHCxp_decision(PPES_rows,3) - mean(ppesxHCxp_decision(PPES_rows,3));
        
        %        oldpes(1:length(oldpes),3) = oldpes(1:length(oldpes),3) - mean(oldpes(1:length(oldpes),3));
        pes(:,3) = pes(:,3) - mean(pes(:,3));
        ppes(PPES_rows,3) = ppes(PPES_rows,3) - mean(ppes(PPES_rows,3));
        %        ppes2(:,3) = ppes2(:,3) - mean(ppes2(:,3));
        calcpes(1:length(calcpes),3) = calcpes(1:length(calcpes),3) - mean(calcpes(1:length(calcpes),3));
        values(1:length(values),3) = values(1:length(values),3) - mean(values(1:length(values),3));
        values_abs(1:length(values_abs),3) = values_abs(1:length(values_abs),3) - mean(values_abs(1:length(values_abs),3));
        
        values_dec = values_dec - mean(values_dec);
        values_dec_abs = values_dec_abs - mean(values_dec_abs);
        
        
        pes_decision = pes_decision - mean(pes_decision([2:48 50:96 98:144 146:192]));
        pes_decxBA = pes_decxBA - mean(pes_decxBA([2:48 50:96 98:144 146:192]));
        pes_decxGB = pes_decxGB - mean(pes_decxGB([2:48 50:96 98:144 146:192]));
        pes_decxG0 = pes_decxG0 - mean(pes_decxG0([2:48 50:96 98:144 146:192]));
        pes_decxHC = pes_decxHC - mean(pes_decxHC([2:48 50:96 98:144 146:192]));
        pes_decxpdecision = pes_decxpdecision - mean(pes_decxpdecision([2:48 50:96 98:144 146:192]));
        pes_decxblockOrder = pes_decxblockOrder - mean(pes_decxblockOrder([2:48 50:96 98:144 146:192]));
        
        pes_decxBAxpdecision = pes_decxBAxpdecision - mean(pes_decxBAxpdecision([2:48 50:96 98:144 146:192]));
        pes_decxGBxpdecision = pes_decxGBxpdecision - mean(pes_decxGBxpdecision([2:48 50:96 98:144 146:192]));
        pes_decxG0xpdecision = pes_decxG0xpdecision - mean(pes_decxG0xpdecision([2:48 50:96 98:144 146:192]));
        pes_decxHCxpdecision = pes_decxHCxpdecision - mean(pes_decxHCxpdecision([2:48 50:96 98:144 146:192]));
        pes_decxBAxOrder = pes_decxBAxOrder - mean(pes_decxBAxOrder([2:48 50:96 98:144 146:192]));
        
        ppes_decision = ppes_decision - mean(ppes_decision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxBA = ppes_decxBA - mean(ppes_decxBA(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxGB = ppes_decxGB - mean(ppes_decxGB(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxG0 = ppes_decxG0 - mean(ppes_decxG0(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxHC = ppes_decxHC - mean(ppes_decxHC(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxpdecision = ppes_decxpdecision - mean(pes_decxpdecision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxBAxpdecision = ppes_decxBAxpdecision - mean(ppes_decxBAxpdecision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxGBxpdecision = ppes_decxGBxpdecision - mean(ppes_decxGBxpdecision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxG0xpdecision = ppes_decxG0xpdecision - mean(ppes_decxG0xpdecision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        ppes_decxHCxpdecision = ppes_decxHCxpdecision - mean(ppes_decxHCxpdecision(intersect([2:48 50:96 98:144 146:192],find(ppes_rows_dec==1))));
        
        %        figure(96); clf;
        %        subplot(3,1,1);
        %        hold on;
        %       plot(oldpes(1:48,3),'b');                             %OLD PES
        %       plot(old_old(1:48,3),'r');
        %        plot(pes(1:48,3),'k');                                 %MODEL-TRACKED PES
        %        subplot(3,1,2);
        %        plot(values(1:48,3),'r');                             %MODEL-tracked VALUE
        %plot(circshift(values(1:48), total_trials+1),'r');   %VALUE ALIGNED W/ DECISIONS MADE
        %        subplot(3,1,3);
        %        hold on;
        %        plot(subject.b.decisions(1:48), 'b');
        %        plot(trustee_action(1:48),'k');
        %        hold off;
        %
        
        
        %% writing PEs
        %           dlmwrite([write_location '/trust' id 'values' '.dat'], values,'delimiter','\t','precision','%.6f');
        %           dlmwrite([write_location '/trust' id 'values_abs' '.dat'], values_abs,'delimiter','\t','precision','%.6f');
        
        %           dlmwrite([write_location '/trust' id 'values_dec' '.dat'],[regressors(:, 1:2) values_dec],'delimiter','\t','precision','%.6f');
        %           dlmwrite([write_location '/trust' id 'values_abs_dec' '.dat'],[regressors(:, 1:2) values_dec_abs],'delimiter','\t','precision','%.6f');
        
        %            dlmwrite([write_location '/trust' id 'PPEs' '.dat'], ppes,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'BAxPPEs' '.dat'],ppesxBA,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'GBxPPEs' '.dat'],ppesxGB,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'G0xPPEs' '.dat'],ppesxG0,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'HCxPPEs' '.dat'],ppesxHC,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxp_decision' '.dat'],ppesxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxBAxp_decision' '.dat'],ppesxBAxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxGBxp_decision' '.dat'],ppesxGBxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxG0xp_decision' '.dat'],ppesxG0xp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxHCxp_decision' '.dat'],ppesxHCxp_decision,'delimiter','\t','precision','%.6f');
        
        %            dlmwrite([write_location '/trust' id 'PPEs_d' '.dat'], [regressors(:, 1:2) ppes_decision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'BAxPPEs_d' '.dat'],[regressors(:, 1:2) ppes_decxBA],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'GBxPPEs_d' '.dat'],[regressors(:, 1:2) ppes_decxGB],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'G0xPPEs_d' '.dat'],[regressors(:, 1:2) ppes_decxG0],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'HCxPPEs_d' '.dat'],[regressors(:, 1:2) ppes_decxHC],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxp_decision_d' '.dat'],[regressors(:, 1:2) ppes_decxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxBAxp_decision_d' '.dat'],[regressors(:, 1:2) ppes_decxBAxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxGBxp_decision_d' '.dat'],[regressors(:, 1:2) ppes_decxGBxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxG0xp_decision_d' '.dat'],[regressors(:, 1:2) ppes_decxG0xpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PPEsxHCxp_decision_d' '.dat'],[regressors(:, 1:2) ppes_decxHCxpdecision],'delimiter','\t','precision','%.6f');
        %
        dlmwrite([write_location '/trust' id 'PEs' '.dat'], pes,'delimiter','\t','precision','%.6f');
        %             dlmwrite([write_location '/trust' id 'BAxPEs' '.dat'],pesxBA,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'GBxPEs' '.dat'],pesxGB,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'G0xPEs' '.dat'],pesxG0,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'HCxPEs' '.dat'],pesxHC,'delimiter','\t','precision','%.6f');
        %             dlmwrite([write_location '/trust' id 'PEsxtrial' '.dat'],[regressors(:,1:2) pesxTrial],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxp_decision' '.dat'],pesxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxBAxp_decision' '.dat'],pesxBAxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxGBxp_decision' '.dat'],pesxGBxp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxG0xp_decision' '.dat'],pesxG0xp_decision,'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxHCxp_decision' '.dat'],pesxHCxp_decision,'delimiter','\t','precision','%.6f');
        
        %           dlmwrite([write_location '/trust' id 'PEs_d' '.dat'],[regressors(:, 1:2) pes_decision],'delimiter','\t','precision','%.6f');
        %           dlmwrite([write_location '/trust' id 'BAxPEs_d' '.dat'],[regressors(:, 1:2) pes_decxBA],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'GBxPEs_d' '.dat'],[regressors(:, 1:2) pes_decxGB],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'G0xPEs_d' '.dat'],[regressors(:, 1:2) pes_decxG0],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'HCxPEs_d' '.dat'],[regressors(:, 1:2) pes_decxHC],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxpdecision_d' '.dat'],[regressors(:, 1:2) pes_decxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxBAxpdecision_d' '.dat'],[regressors(:, 1:2) pes_decxBAxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxGBxpdecision_d' '.dat'],[regressors(:, 1:2) pes_decxGBxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxG0xpdecision_d' '.dat'],[regressors(:, 1:2) pes_decxG0xpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxHCxpdecision_d' '.dat'],[regressors(:, 1:2) pes_decxHCxpdecision],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'PEsxOrder_d' '.dat'],[regressors(:, 1:2) pes_decxblockOrder],'delimiter','\t','precision','%.6f');
        %            dlmwrite([write_location '/trust' id 'BAxPEsxOrder_d' '.dat'],[regressors(:, 1:2) pes_decxBAxOrder],'delimiter','\t','precision','%.6f');
        
    end
    catch
        disp([filename ' something went wrong'])
        continue
    end
    close all;
end

%collinearity
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