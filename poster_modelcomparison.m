
%simplest model
L_c0k0m0r0h0p0n0a0u0 = load('L_counter0_multisession0_fixed1_SigmaKappa0_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret0');

%model with kappa
L_c0k1m0r0h0p0n0a0u0 = load('L_counter0_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret0');

Lgood = [L_c0k0m0r0h0p0n0a0u0.L;L_c0k1m0r0h0p0n0a0u0.L];
Lgood(:,[1:20,22:44]);

[posterior,out] = VBA_groupBMC(Lgood);

%model with kappa + counterfactual
L_c1k1m0r0h0p0n0a0u0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret0');

Lgood = [L_c0k1m0r0h0p0n0a0u0.L; L_c1k1m0r0h0p0n0a0u0.L];
Lgood(:,[1:20,22:44]); %poster subjects
Lgood = Lgood(:,[1:23,25:52,54]); %all subjects
[posterior,out] = VBA_groupBMC(Lgood);



%model with kappa + multisession
L_c1k1m1r0h0p0n0a0u0 = load('L_counter1_multisession1_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret0');
Lgood = [L_c1k1m0r0h0p0n0a0u0.L; L_c1k1m1r0h0p0n0a0u0.L];
Lgood=Lgood(:,[1:20,22:44]); %poster subjects
[posterior,out] = VBA_groupBMC(Lgood);


