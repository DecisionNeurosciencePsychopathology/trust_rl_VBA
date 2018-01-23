function param_recov_struct = run_trust_parameter_recovery
%Script is a wrapper to run the parameter recovery script over a model
%group

%Get all the subjects from a model specific directory
%subdirs = {'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT' 'new_f_trust_Qlearn_counter_hybrid_regret_pmv' 'new_f_trust_Qlearn_counter_trustee' 'new_f_trust_Qlearn_null_pmv'};
models = {'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT', 'f_trust_Qlearn_counter_hybrid'};
models={'f_trust_Qlearn_counter_hybrid'};

%Set up data ouput struct
param_recov_struct = struct();


%TODO - pull and save all the unique f and g names, save those as evo and
%oversvation pairs for model confusion, then re-run this script with the
%caveat that is model confusion is on substitue f and g functions, need to
%think of the output though...

%Look up table idea
%If lookup table D.N.E. create the empty table
%cols are model name f_fname g_fname
%if table D.E. check if model is already in lookup table
%to run model confusion table height must be > 2

for k = 1:length(models)
    subj_mat_files = glob([models{k} '/*.mat']);
    for j = 1:length(subj_mat_files) %Loop trhough the subjects
        load(subj_mat_files{j}) %i got saved for whatever reason
        close all; %Remove graphics if any appear
        try
            [param_recov_struct.(models{k}).post{j,1}, param_recov_struct.(models{k}).out{j,1}] = trust_param_recovery(out,j);
        catch
            fprintf('\nVBA broke down probably due to not having all the proper functions\n')
        end
    end
end



function [post_new, out_new] = trust_param_recovery(out,seed)
%Create psudorandom datasets from the VBA model output
%We need the 'out' from VBA_NLStateSpaceModel, specifically
%out.stuffStat.gx
%seed  = randomseed for creating pseudo data points

%TODO 0 add fname and gname for model confusion step

%Set up seed
rng(seed)

%Start par loop if not set up already
p = gcp('nocreate');
if isempty(p)
    parpool('local')
end

%Pull the required info from the out struct
f_fname = out.options.f_fname; %Evolution function (Learning rule)
g_fname = out.options.g_fname; %Observation function (Choice Rule)
u = out.u; %Input info (rew stuct, ect)
gx = out.suffStat.gx; %The simulated data
dim = out.dim; %Dimensions of hiddens states and parameters
options = out.options; %All previosuly defined options

%Create n pseudo subjects based on one subejcts model responses thorugh
%random sampling
num_samples = 3;
y = gen_psudo_subj_responses(num_samples,gx);
parfor i = 1:num_samples
    [post_new(i),out_new(i)] = VBA_NLStateSpaceModel(y(i,:),u,f_fname,g_fname,dim,options);
end


function y = gen_psudo_subj_responses(num_samples,model_responses)

num_resps = length(model_responses);
for i = 1:num_samples
    y(i,:) = double(rand(num_resps,1)' < model_responses);
end

