%Script will extract and right all PEs for given models (i.e. subdirs)

%Which models to use
%subdirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_corrected', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};
subdirs={'f_trust_Qlearn1'};

%Extract the PEs
trust_extract_PEs(subdirs)

%Really, really need to fix this cd stuff, it really is a pain to deal
%with!
current_dir = pwd; %This should be trust_rl_VBA

%Write the PEs
for subdir = subdirs
    cd(current_dir)
    trust_writeregs_PEs(subdir)
    cd(current_dir)
    asterisk(subdir) %Run asterisk script to handle the blocks
end


