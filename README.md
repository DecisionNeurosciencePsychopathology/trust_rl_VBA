# trust_rl_VBA
VBA-style Q-learning and other models of modified Trust Game
1. trust_fit_group_vba.m: gets data from all subjects (.mat), runs a Qlearning (trust_Qlearning_ushifted) algorithm w/ preset parameters, saves F(it)-values in table.
2. trust_Qlearning_ushifted.m: VBA fitting of Qlearning to a single subject trust data. Note that the u vector (which includes the subject's behavior +outcomes+experimental design) is shifted by one trial ahead, because enivronmental feedback on trial t-1 is what determines the value of the hidden states on trial t. Calls to f() and g() (obervation/evolution functions, respectively). Saves the figure and the results of the VBA inversion routine.
3. f_trust_Qlearn_counter.m: evolution function where reinforcement reflects counterfactual rewards (with respect to the share action)
4. f_trust_Qlearn1.m: simpler evolution function where reinforcement reflects counterfactual feedback from trustee (also with respect 
to the share actions).
5. f_trust_Qlearn2.m: simpler evolution function where reinforcement for keep/share is considered separately.
6. g_trust_softmax_ED: observation function that includes the kappa parameter (bias for keep/share actions).
7. g_trust_softmax_1: simpler observation function.
8. trust_extract_PEs: gets the PEs from saved .mat files resulting from VBA inversion routines, saves in table.
9. trust_write_PEs: writes regressors based on the table resulting from trust_extract_PEs.
10. trust_writevalues: gets and writes value regressors.
11. trust_modelComparisons: helper code to compares F(its) between models. Uses saved tables from trust_fit_group_vba.m.
