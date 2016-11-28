%1-Uncorrected
%2-Corrected
%3-Trustee
%4-Non-Counterfactual
%5-Hybrid
%6-Hybrid Regret
%7-Trustee Hybrid

L_behav=[];
L_scan=[];
% load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\filtered\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
% L_scan(1,:) = L;
% load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\corrected\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
% L_scan(2,:) = L;
load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\trustee\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
L_scan(3,:) = L;
%load('E:\trust_model_comparision\trust_rl_VBA\all_behav\all_ages\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
% L_behav(1,:)=L;
% load('E:\trust_model_comparision\trust_rl_VBA\all_behav\corrected\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
% L_behav(2,:)=L;
load('E:\trust_model_comparision\trust_rl_VBA\all_behav\trustee\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
L_behav(3,:)=L;

load('E:\trust_model_comparision\trust_rl_VBA\all_behav\all_ages\L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
L_behav(4,:)=L;
load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\filtered\L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
L_scan(4,:) = L;

load('E:\trust_model_comparision\trust_rl_VBA\all_behav\hybrid\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
L_behav(5,:)=L;
load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\hybrid\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
L_scan(5,:) = L;

load('E:\trust_model_comparision\trust_rl_VBA\all_behav\regret\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
L_behav(6,:)=L;
load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\regret\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
L_scan(6,:) = L;


% load('E:\trust_model_comparision\trust_rl_VBA\all_behav\trustee_hybrid\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
% L_behav(7,:)=L;
% load('E:\trust_model_comparision\trust_rl_VBA\hc_scan\trustee_hybrid\L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
% L_scan(7,:) = L;

[post_scan,out_scan] = VBA_groupBMC(L_scan(3:6,:))
[post_behav,out_behav] = VBA_groupBMC(L_behav(3:6,:))