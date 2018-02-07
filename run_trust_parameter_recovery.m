function [param_recov_struct, model_confusion_struct] = run_trust_parameter_recovery
%Script is a wrapper to run the parameter recovery script over a model
%group

%Get all the subjects from a model specific directory
%subdirs = {'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT' 'new_f_trust_Qlearn_counter_hybrid_regret_pmv' 'new_f_trust_Qlearn_counter_trustee' 'new_f_trust_Qlearn_null_pmv'};
models = {'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT', 'f_trust_Qlearn_counter_hybrid'};
models={'f_trust_Qlearn_counter_hybrid'};
models = {'actual','f_trust_SVM1','policy','regret','trustee'};

%Set up path to data
data_path = 'E:/Box Sync/skinner/projects_analyses/Project Trust/data/model-derived/scan/';

%Set up data ouput struct
param_recov_struct = struct();

%Model confusion lookup table
model_confusion_flag = 1; %for now hard code

%Create/load model confusion lookup table
if ~exist('model_lookup.mat', 'file') && model_confusion_flag
    model_lookup = create_model_lookup_table(models,data_path);
    save('model_lookup','model_lookup')
else
    load('model_lookup') %var name is "model_lookup"
end


%Look up table idea
%If lookup table D.N.E. create the empty table
%cols are model name f_fname g_fname
%if table D.E. check if model is already in lookup table
%to run model confusion table height must be > 2

tic
for k = 1:length(models)
    subj_mat_files = glob([data_path models{k} '/*.mat']);
    for j = 1:length(subj_mat_files) %Loop through the subjects
        load(subj_mat_files{j}) %i got saved for whatever reason
        close all; %Remove graphics if any appear
        try
            f_fname = out.options.f_fname; %Evolution function (Learning rule)
            g_fname = out.options.g_fname; %Observation function (Choice Rule)
            [param_recov_struct.(models{k}).post{j,1}, param_recov_struct.(models{k}).out{j,1}, pseudo_y] = trust_param_recovery(out, dim, options, f_fname, g_fname, j);
            param_recov_struct.(models{k}).psudo_y{j,1} = pseudo_y; %Save the pseudo y
            
            if model_confusion_flag
                models_to_run = find(~strcmp(models{k},model_lookup.model_name));
                for m = models_to_run'
                    [model_confusion_struct.([models{k} '_confused_with_' model_lookup.model_name{m}]).post{j,1}, model_confusion_struct.([models{k} '_confused_with_' model_lookup.model_name{m}]).out{j,1}] = ...
                        trust_param_recovery(out, model_lookup.dim{m}, model_lookup.options{m}, model_lookup.f_fname{m}, model_lookup.g_fname{m}, j, pseudo_y);
                    model_confusion_struct.([models{k} '_confused_with_' model_lookup.model_name{m}]).pseudo_y{j,1} = pseudo_y; %Duplicate but just as a sanity check
                end
            end
        catch
            fprintf('\nVBA broke down probably due to not having all the proper functions\n')
        end
    end
end
toc


function [post_new, out_new, y] = trust_param_recovery(out, dim, options, f_fname, g_fname, seed, y)
%Create psudorandom datasets from the VBA model output
%We need the 'out' from VBA_NLStateSpaceModel, specifically
%out.stuffStat.gx
%seed  = randomseed for creating pseudo data points

try y; catch, y=[]; end;

%Set up seed
rng(seed)

%Turn off figure creation
options.DisplayWin = 0;

%Start par loop if not set up already
p = gcp('nocreate');
if isempty(p)
    parpool('local')
end

%Pull the required info from the out struct
%f_fname = out.options.f_fname; %Evolution function (Learning rule)
%g_fname = out.options.g_fname; %Observation function (Choice Rule)
u = out.u; %Input info (rew stuct, ect)
gx = out.suffStat.gx; %The simulated data
%dim = out.dim; %Dimensions of hiddens states and parameters
%options = out.options; %All previosuly defined options

%Create n pseudo subjects based on one subejcts model responses thorugh
%random sampling, if y is not supplied initially
num_samples = 3;
if isempty(y)
    y = gen_psudo_subj_responses(num_samples,gx);
end

parfor i = 1:num_samples
    [post_new(i),out_new(i)] = VBA_NLStateSpaceModel(y(i,:),u,f_fname,g_fname,dim,options);
end


function y = gen_psudo_subj_responses(num_samples,model_responses)
%Simulate subject model responses using a weighted binomial distribution

%Initialize y
y = zeros(num_samples,length(model_responses));

%Used weighted binomial sampling.
for i = 1:num_samples
    y(i,:) = binornd(1,model_responses); 
end

function model_lookup = create_model_lookup_table(models,data_path)

%Initialize tables
model_lookup = table();
tmp_table = table();

for model = models
    %Just load in one subject file per model
    subj_mat_files = glob([data_path model{:} '/*.mat']);
    load(subj_mat_files{1});
    close all; %Pesky figures
    
    tmp_table.model_name = model;
    tmp_table.f_fname = {f_fname};
    tmp_table.g_fname = {g_fname};
    tmp_table.dim = {dim};
    tmp_table.options = {options};
    model_lookup = [model_lookup; tmp_table];
end
