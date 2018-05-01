function bmc_output=run_model_confusion_analysis
%Script will load in the model confusion matrix and reshape the log model
%evidence to make one "L array" per model.

%Load in the model confusion var (At this point we still in the original
%parameter recovery simulations)
load('E:\trust_model_comparision\trust_rl_VBA\model_confusion/model_conf_output.mat') %model_confusion_struct
load('E:\trust_model_comparision\trust_rl_VBA\trust_param_recov_output.mat') %param_recov_struct

%Define models
models = {'actual','f_trust_SVM1','policy','regret','trustee'};
bmc_output = struct();

%Pull the data apart
data=compile_model_confusion_data(model_confusion_struct,param_recov_struct,models);

%GroupBMC
for i = 1:length(models)    
    [bmc_output.(models{i}).posterior, bmc_output.(models{i}).out,h] = run_model_confusion_groupBMC(data{i});
    bmc_output.(models{i}).ep = bmc_output.(models{i}).out.ep; %Estimated probabilites
    bmc_output.(models{i}).ef = bmc_output.(models{i}).out.Ef; %Estimated frequencies
    savefig(h,[models{i} 'groupBMC']) %Save figure
end


%% Sub Functions
function data=compile_model_confusion_data(model_confusion_struct,param_recov_struct,models)
fnames = fieldnames(model_confusion_struct);
for fname = fnames'
    original_model = regexp(fname{:},'(\w+)_c','tokens');
    model_index=~cellfun(@isempty,strfind(models,original_model{:}));
    [x,y]=size(model_confusion_struct.(fname{:}).L);
    n_col = x*y;
     if ~exist('data')
         data = cell(length(models),1);
     end
    data{model_index} = [data{model_index,:}; reshape(model_confusion_struct.(fname{:}).L,1,n_col)];
end

%Add in the parameter recovery log model evidence
for i = 1:length(models)
    [x,y]=size(param_recov_struct.(models{i}).L);
    n_col = x*y;
    data{i} = [data{i,:}; reshape(param_recov_struct.(models{i}).L,1,n_col)];
end


%% ------------------------------------------------------------------------
function [posterior,out,h] = run_model_confusion_groupBMC(data)
close all; %Remove the graphics
[posterior,out] = VBA_groupBMC(data);
h = figure(1);