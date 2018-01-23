%Script will extract and right all PEs for given models (i.e. subdirs)

%Which models to use
%subdirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_corrected', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};
subdirs={'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT','new_f_trust_Qlearn_counter_hybrid_regret_pmv','new_f_trust_Qlearn_counter_trustee','new_f_trust_Qlearn_null_pmv'};
subdirs={'new_f_trust_Qlearn_null_pmv'};
subdirs={'f_trust_Qlearn_policy_censor0_mltrun1_kS_kT'};

%Extract the PEs
trust_extract_PEs(subdirs)

%Really, really need to fix this cd stuff, it really is a pain to deal
%with!
current_dir = pwd; %This should be trust_rl_VBA

%Write the PEs
for subdir = subdirs    
    cd(current_dir)
    just_extract_model_pes_condensed(subdir)
    %trust_writeregs_PEs(subdir)
    cd(current_dir)
    asterisk(subdir) %Run asterisk script to handle the blocks
end


