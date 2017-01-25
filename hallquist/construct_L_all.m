load('E:\trust_model_comparision\trust_rl_VBA\hallquist\L_null_hall_n_30.mat')
l_hall_struct.null = L;
load('E:\trust_model_comparision\trust_rl_VBA\hallquist\L_hybrid_hall_n_30.mat')
l_hall_struct.hybrid = L;
load('E:\trust_model_comparision\trust_rl_VBA\hallquist\L_regret_hall_n_30.mat')
l_hall_struct.hybrid_regret = L;
load('E:\trust_model_comparision\trust_rl_VBA\hallquist\L_trustee_hall_n_30.mat')
l_hall_struct.trustee = L;

l_hall_struct.L_all(1,:) = l_hall_struct.null;
l_hall_struct.L_all(2,:) = l_hall_struct.hybrid;
l_hall_struct.L_all(3,:) = l_hall_struct.hybrid_regret;
l_hall_struct.L_all(4,:) = l_hall_struct.trustee;

[posterior,out] = VBA_groupBMC(l_hall_struct.L_all);


%Linas behav subjects for scan use compare_all_models.m
load('E:\trust_model_comparision\trust_rl_VBA\final_behavior\L_null_behav_n_15.mat')
l_behav_struct.null = L;
load('E:\trust_model_comparision\trust_rl_VBA\final_behavior\L_hybrid_behav_n_15.mat')
l_behav_struct.hybrid = L;
load('E:\trust_model_comparision\trust_rl_VBA\final_behavior\L_hybrid_regret_behav_n_15.mat')
l_behav_struct.hybrid_regret = L;
load('E:\trust_model_comparision\trust_rl_VBA\final_behavior\L_trustee_behav_n_15.mat')
l_behav_struct.trustee = L;
load('E:\trust_model_comparision\trust_rl_VBA\final_behavior\L_subjCounter_behav_n_15.mat')
l_behav_struct.subjCounter = L;

l_behav_struct.L_all(1,:) = l_behav_struct.null;
l_behav_struct.L_all(2,:) = l_behav_struct.hybrid;
l_behav_struct.L_all(3,:) = l_behav_struct.hybrid_regret;
l_behav_struct.L_all(4,:) = l_behav_struct.trustee;
l_behav_struct.L_all(5,:) = l_behav_struct.subjCounter;

[posterior,out] = VBA_groupBMC(l_behav_struct.L_all);