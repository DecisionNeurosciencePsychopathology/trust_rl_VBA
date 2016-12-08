# trust_rl_VBA
VBA-style Q-learning and other models of modified Trust Game.
A procedure to run programs in this collection can go as follows. (1) run trust_fit_group_vba. (2) After determining desired model, run trust_extract_PEs. (3) Run trust_write_PEs and trust_writevalues.

1. trust_fit_group_vba.m: gets data from all subjects (.mat), runs a Qlearning (trust_Qlearning_ushifted) algorithm w/ preset parameters, saves F(it)-values in table.

2a. trust_Qlearning_ushifted.m: VBA fitting of Qlearning to a single subject trust data. Note that the u vector (which includes the subject's behavior +outcomes+experimental design) is shifted by one trial ahead, because enivronmental feedback on trial t-1 is what determines the value of the hidden states on trial t. Calls to f() and g() (obervation/evolution functions, respectively). Saves the figure and the results of the VBA inversion routine.

2b. trust_Qlearning_ushifted2.m: implemented for using with Social Value Model aka Fareri et al., 2015.

3. Evolution functions:
  3a. f_trust_Qlearn_counter.m: reinforcement reflects counterfactual rewards (with respect to the share action).
  3b. f_trust_Qlearn1.m: reinforcement reflects counterfactual feedback from trustee (also with respect to the share actions).
  3c. f_trust_Qlearn2.m: simpler evolution function where reinforcement for keep/share is considered separately.
  3d. f_trust_Qlearn_counter_corrected.m: subject-counterfactual = actual rewards - would-be outcome from subject's alternative action    
  3e. f_trust_Qlearn_counter_hybrid.m: subject-counterfactual *counterfactuals are only calculated when participant's action = KEEP*.
  3e. f_trust_Qlearn_counter_trustee.m: trustee-counterfactual = actual rewards - would-be outcome from trustee's alternative action. 
  3f. f_trust_Qlearn_counter_hybrid_regret.m: regret = actual rewards - maximum outcome (1.5).
  3g. f_trust_Qlearn_counter_hybrid_regret_purist: regret *only when there is something to regret*
  3h. f_trust_SVM1: Social Value Model aka Fareri et al., 2015.

4. Observation functions:
  4a. g_trust_softmax_ED: observation function that includes the kappa parameter (bias for keep/share actions).
  4b. g_trust_softmax_1: simpler observation function.
  4c. g_trust_SVM1.m: Social Value Model aka Fareri et al., 2015 *also softmax*.

5. trust_extract_PEs: gets the PEs from saved .mat files resulting from VBA inversion routines, saves in table.

6. trust_write_PEs: writes regressors based on the table resulting from trust_extract_PEs.

7. trust_writevalues: gets and writes value regressors.

8. trust_modelComparisons: helper code to compares F(its) between models. Uses saved tables from trust_fit_group_vba.m.
