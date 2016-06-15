L_p0n1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1.mat');
L_p1n0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0.mat');
L_p1n1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n1.mat');

L_m0r0h0p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h0p0n1a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1_assymetry_choice0.mat');
L_m0r0h0p1n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0_assymetry_choice0.mat');
L_m0r1h0p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation1_humanity0_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h1p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity1_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h0p1n1a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n1_assymetry_choice0.mat');
L_m0r0h0p0n0a1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1.mat');
L_k0m0r0h0p0n0a1 = load('L_multisession0_fixed1_SigmaKappa0_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1.mat');



Lgood = [L_m0r0h0p0n0a0.L;L_m0r1h0p0n0a0.L;L_m0r0h1p0n0a0.L;L_m0r0h0p0n1a0.L;L_m0r0h0p1n0a0.L;L_m0r0h0p1n1a0.L; L_m0r0h0p0n0a1.L; L_k0m0r0h0p0n0a1.L];


%%
clear all; close all;
new = load('Ushifted2_L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lnew = new.L;
old = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lold = old.L;

[posterior,out] = VBA_groupBMC([Lnew; Lold]);

clear all; close all;
old = load('L_counter0_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lold = old.L;

new = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lnew = new.L;
Lgood = [Lnew; Lold];
Lgood = Lgood(:,[1:10,12:28,30:47,49]);

[posterior,out] = VBA_groupBMC([Lnew; Lold]);

L_c1m0r0h0p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m1r0h0p0n0a0 = load('L_counter1_multisession1_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r1h0p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation1_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h1p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity1_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h0p1n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h0p0n1a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1_assymetry_choice0_beta0');
L_c1m0r0h0p0n0a1 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1_beta0');

Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L;L_c1m0r0h0p1n0a0.L;L_c1m0r0h0p0n1a0.L;L_c1m0r0h0p0n0a1.L];

[posterior,out] = VBA_groupBMC(Lgood);

Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L];

[posterior,out] = VBA_groupBMC(Lgood);

Lgood = [L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L];

[posterior,out] = VBA_groupBMC(Lgood([1,2],:));

Lgood = [L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L;L_c1m0r0h0p1n0a0.L;L_c1m0r0h0p0n1a0.L;L_c1m0r0h0p0n0a1.L];

L_c1m0r0h0p0n0a0u1 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret1');
Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r0h0p0n0a0u1.L];
[posterior,out] = VBA_groupBMC(Lgood(:,[1:10,12:28,30:47,49])); %L_c1m0r0h0p0n0a0.L > L_c1m0r0h0p0n0a0u1.L