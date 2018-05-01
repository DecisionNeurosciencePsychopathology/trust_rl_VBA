%Real quick way to make the difference between hybrid model PE's and null
%model PE's, go into the trust_wrtieregs_PEs.m script and replace the
%typical M_PE's file name with the new difference name, then rewrite the
%script to save it as the diff PE regressor.

load('E:\trust_model_comparision\trust_rl_VBA\f_trust_Qlearn1\PEs\modelPEs_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
null_m = M;
load('E:\trust_model_comparision\trust_rl_VBA\f_trust_Qlearn_counter_hybrid\PEs\modelPEs_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat')
hybrid_m = M;
clear M
M(:,1) = null_m(:,1);
M(:,2:193) = hybrid_m(:,2:end) - null_m(:,2:end);