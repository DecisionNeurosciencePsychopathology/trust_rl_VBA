function recov_output=run_parameter_recovery_analysis
%Load in the paramter recovery file
load('trust_param_recov_output.mat') %param_recov_struct

%Per model per parameter
%x: n_subs X 100 clones
%y: n_subs X 100 psudo samples

%Grab model names
model_names = fieldnames(param_recov_struct);

%Do we want to save the plotting data for R?
save_for_plotting = 1;

%Initialize final output struct
recov_output = struct();

%Point to where the single subject data is
original_data_path = 'E:/Box Sync/skinner/projects_analyses/Project Trust/data/model-derived/scan/';

        ct = 1;

for model = model_names' %Model loop
    
    %For saving later
    T = table();
    
    %Compile the original data
    parameter_table=compile_vba_parameters([original_data_path model{:}]);
    
    %Create the x axis variable
    %Loop through all vaiables in paramter table?
    original_params = struct();
    params = fieldnames(parameter_table);
    params(strcmp(params(:,1), 'ID'),:) = []; %Remove ID & Properties
    params(strcmp(params(:,1), 'Properties'),:) = [];
    
    for param = params'
        original_params.(param{:}) = repmat(parameter_table.(param{:}),1,100);
    end
    
    %Create the y axis variable
    simulated_params = struct();    
    simulated_params=return_compiled_simulated_params(param_recov_struct,model,'theta',simulated_params);
    simulated_params=return_compiled_simulated_params(param_recov_struct,model,'phi',simulated_params);
    
    
    %We have to rename the phi parameters from 3-6 because they are not
    %fixed for both the simulated and originial parameter sets.
    if ~strcmpi(model{:},'f_trust_svm1')
        simulated_params = restructure_kappa_parameters(simulated_params, [original_data_path model{:}], model{:});
        original_params = restructure_kappa_parameters(original_params, [original_data_path model{:}], model{:});
        
    else
        simulated_params=rename_parameters(simulated_params, model{:}); 
        original_params=rename_parameters(original_params, model{:}); 
    end  
    
    %Loop thorugh all params to make correlations and plots
    %Could have a lookup table to rename these guys
    vba_params = fieldnames(original_params);
    
    %TODO Check if original and simulated params do not have the same fieldnames
    %and kick out

    %Run the correlations and svae the output per parameter per model
    for vba_param = vba_params'
        if iscell(simulated_params.(vba_param{:}))
            simulated_params.(vba_param{:})=cell2mat(simulated_params.(vba_param{:}));
        end
        [rows,cols] = size(original_params.(vba_param{:}));
        n_cols = rows*cols;
        [recov_output.(model{:}).(vba_param{:}).r,recov_output.(model{:}).(vba_param{:}).p] = ...
            corrcoef(reshape(original_params.(vba_param{:}),n_cols,1),reshape(simulated_params.(vba_param{:}),n_cols,1),'rows','complete');
        
        figure(ct)
        scatter(reshape(original_params.(vba_param{:}),n_cols,1),reshape(simulated_params.(vba_param{:}),n_cols,1))
        title(sprintf('MODEL %s : PARAMETER: %s',model{:},vba_param{:}))
        ct = ct+1;
        
        
        if save_for_plotting
            R_dir = 'for_R_plotting';
            mkdir(R_dir)
            T.(['original_' vba_param{:}]) = reshape(original_params.(vba_param{:}),n_cols,1);
            T.(['simulated_' vba_param{:}]) = reshape(simulated_params.(vba_param{:}),n_cols,1);            
        end
    end
    
    if ~isempty(T) && save_for_plotting
        writetable(T, [R_dir filesep [model{:} 'parameters_for_plotting.csv']])
    end
    
end


function simulated_params=return_compiled_simulated_params(param_recov_struct,model,param_of_interest,simulated_params)


switch param_of_interest
    case 'theta'
        p_str = 'theta_param_';
        n_vars = length(param_recov_struct.(model{:}).muTheta{1,1});
        mu_str = {'muTheta'};
    case 'phi'
        p_str = 'phi_param_';
        n_vars = length(param_recov_struct.(model{:}).muPhi{1,1});
        mu_str = {'muPhi'};
    otherwise
        error('Options are with theta or phi')
end

    for i = 1:n_vars
        tmp_str = {[p_str num2str(i)]};
        [rows,cols] = size(param_recov_struct.(model{:}).(mu_str{:}));
        n_cols = rows*cols;
        tmp=reshape(param_recov_struct.(model{:}).(mu_str{:}).',1,n_cols);
        
        param_lenth = cellfun(@length, tmp);
        max_params = max(cellfun(@length, tmp)); %This assumes that the higest number of parameters is correct!
        
        %Find where there could be missing parameters for subjects and fill them with nans
        missing_param_idx=find(param_lenth<max_params);
        
        %Add in nans is need be
        if ~isempty(missing_param_idx)
            for idx = missing_param_idx
                filler = max_params-length(tmp{idx});
                tmp{idx} = [tmp{idx}; nan(filler,1)];
            end
        end
        simulated_params.(tmp_str{:}) = permute(reshape(cellfun(@(v)v(i),tmp),cols,rows),[2,1]);
    end

    function s = restructure_kappa_parameters(s,path_to_subject_data,model_name)
        %Becasue the trustee kappa parameter is not fixed in the observation
        %function during multisession we have to rename/ restructe our phi
        %matrices
        
        %Grab all the subject files
        subj_files = glob([path_to_subject_data filesep '*.mat']);
        
        for i = 1:length(subj_files)            
            load(subj_files{i},'b') %only load in b
            
            %Hardcoded, but this seems to be the only consistant variable
            trustee_order = {b.identity{1} b.identity{50} b.identity{97} b.identity{173}};                    
            
            %Reconstruction
            trustees = {'good';'bad';'neutral';'computer'};
            for j = 1:length(trustees)
                trustee_idx = find(strcmp(trustee_order, trustees{j}))+2;
                if ~isempty(trustee_idx)
                    s.(['kappa_' trustees{j}])(i,:) = s.(['phi_param_' num2str(trustee_idx)])(i,:);
                else
                    %If the data is missing replace the row with nans
                    s.(['kappa_' trustees{j}])(i,:) = nan(1,length(s.phi_param_1));
                end
            end
        end
        
        s=rename_parameters(s,model_name);


        function s = rename_parameters(s,model_name)
            if ~strcmpi(model_name,'f_trust_svm1')
                %Renaming vars
                s.('learning_rate') = s.theta_param_1;
                s.('beta') = s.phi_param_1;
                s.('kappa_subject') = s.phi_param_2;
                
                %Remove old fields
                fields_to_remove = matlab.lang.makeUniqueStrings(repmat({'phi_param'},1,7));
                fields_to_remove{1}='theta_param_1';                
                s=rmfield(s,fields_to_remove);
            else
                %Renaming vars
                s.('learning_rate') = s.theta_param_1;
                s.('social_value') = s.theta_param_1;
                s.('beta') = s.phi_param_1;
                
                %Remove
                s=rmfield(s,{'theta_param_1','theta_param_2','phi_param_1'});
            end
            
            

%%%OLD code del later

% % %     n_phi = length(param_recov_struct.(model{:}).muPhi{1,1});
% % %     for i = 1:n_phi
% % %         tmp_str = {['phi_param_' num2str(i)]};
% % %         [rows,cols] = size(param_recov_struct.(model{:}).muPhi);
% % %         n_cols = rows*cols; 
% % %         tmp=reshape(param_recov_struct.(model{:}).muPhi.',1,n_cols);
% % %         
% % %         param_lenth = cellfun(@length, tmp);
% % %         max_params = max(cellfun(@length, tmp)); %This assumes that the higest number of parameters is correct!
% % %         
% % %         %Find where there could be missing parameters for subjects and fill them with nans
% % %         missing_param_idx=find(param_lenth<max_params);
% % %         
% % %         %Add in nans is need be
% % %         if ~isempty(missing_param_idx)
% % %             for idx = missing_param_idx
% % %                 filler = max_params-length(tmp{idx});
% % %                 tmp{idx} = [tmp{idx}; nan(filler,1)];
% % %             end
% % %         end
% % %         simulated_params.(tmp_str{:}) = permute(reshape(cellfun(@(v)v(i),tmp),cols,rows),[2,1]);
% % %     end
% % %     
% % %     
% % %     
% % %         for i = 1:n_theta
% % %         % %         tmp_str = {['theta_param_' num2str(i)]};
% % %         % %         simulated_params.(tmp_str{:}) = [param_recov_struct.(model{:}).muTheta];
% % %         
% % %         tmp_str = {['theta_param_' num2str(i)]};
% % %         [rows,cols] = size(param_recov_struct.(model{:}).muTheta);
% % %         n_cols = rows*cols;
% % %         tmp=reshape(param_recov_struct.(model{:}).muTheta.',1,n_cols);
% % %         
% % %         param_lenth = cellfun(@length, tmp);
% % %         max_params = max(cellfun(@length, tmp)); %This assumes that the higest number of parameters is correct!
% % %         
% % %         %Find where there could be missing parameters for subjects and fill them with nans
% % %         missing_param_idx=find(param_lenth<max_params);
% % %         
% % %         %Add in nans is need be
% % %         if ~isempty(missing_param_idx)
% % %             for idx = missing_param_idx
% % %                 filler = max_params-length(tmp{idx});
% % %                 tmp{idx} = [tmp{idx}; nan(filler,1)];
% % %             end
% % %         end
% % %         simulated_params.(tmp_str{:}) = permute(reshape(cellfun(@(v)v(i),tmp),cols,rows),[2,1]);
% % %     end