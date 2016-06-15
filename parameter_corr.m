array =out2_good;
correlations = zeros(length(array),3);

for subject =1:length(array)
    correlations(subject,:) = [array(subject).diagnostics.C(1,2), array(subject).diagnostics.C(4,1),array(subject).diagnostics.C(4,2)];
end
mean_corr = mean(correlations);
Phi1_Phi2_X0 = zeros(length(out2),3);
for subject =1:length(array)
    Phi1_Phi2_X0(subject,:) = [array(subject).muPhi(1), array(subject).muPhi(2),array(subject).muX0];
end


datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');

%% choose model's parameters
counter = 1;
multisession = 0;
fixed_params_across_runs = 1;
sigma_kappa = 1;
reputation_sensitive = 0;
humanity = 0;
valence_p = 0;
valence_n = 0;
assymetry_choices = 0;
regret = 1;

cd(datalocation{1});
files = dir(strcat('*',sprintf('counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d_regret%d', counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices,regret),'.mat'));
num_of_subjects = length(files);

%omega1
correlations = zeros(num_of_subjects,6);
for ct = 1:num_of_subjects
    filename=files(ct).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    load(filename);
    correlations(ct,:) = [id, out.diagnostics.C(4,1), out.diagnostics.C(4,2), out.diagnostics.C(4,3), out.diagnostics.C(4,5), out.diagnostics.C(4,6)];
%    N(ct,:) = [id, out.suffStat.muX];
end
mean(correlations(:,2))
min(correlations(:,2))
max(correlations(:,2))
histogram(correlations(:,2))

mean(correlations(:,3))
min(correlations(:,3))
max(correlations(:,3))
histogram(correlations(:,3))

%theta(1) = learning rate
mean(correlations(:,4)) %correlated -.35, max = -.67, min = -.009; flat distribution;
min(correlations(:,4))
max(correlations(:,4))
histogram(correlations(:,4))

mean(correlations(:,5))
min(correlations(:,5))
max(correlations(:,5))
histogram(correlations(:,5))

mean(correlations(:,6))
min(correlations(:,6))
max(correlations(:,6))
histogram(correlations(:,6))

%omega2
correlations = zeros(num_of_subjects,6);
for ct = 1:num_of_subjects
    filename=files(ct).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    load(filename);
    correlations(ct,:) = [id, out.diagnostics.C(5,1), out.diagnostics.C(5,2), out.diagnostics.C(5,3), out.diagnostics.C(5,4), out.diagnostics.C(5,6)];
%    N(ct,:) = [id, out.suffStat.muX];
end
mean(correlations(:,2))
min(correlations(:,2))
max(correlations(:,2))
histogram(correlations(:,2))

mean(correlations(:,3))
min(correlations(:,3))
max(correlations(:,3))
histogram(correlations(:,3))

%theta(1) = learning rate
mean(correlations(:,4)) %correlated -.35, max = -.67, min = -.009; flat distribution;
min(correlations(:,4))
max(correlations(:,4))
histogram(correlations(:,4))

mean(correlations(:,5))
min(correlations(:,5))
max(correlations(:,5))
histogram(correlations(:,5))

mean(correlations(:,6))
min(correlations(:,6))
max(correlations(:,6))
histogram(correlations(:,6))

%single omega

correlations = zeros(num_of_subjects,5);
for ct = 1:num_of_subjects
    filename=files(ct).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    load(filename);
    correlations(ct,:) = [id, out.diagnostics.C(4,1), out.diagnostics.C(4,2), out.diagnostics.C(4,3), out.diagnostics.C(4,5)];
%    N(ct,:) = [id, out.suffStat.muX];
end
mean(correlations(:,2))
min(correlations(:,2))
max(correlations(:,2))
histogram(correlations(:,2))

mean(correlations(:,3))
min(correlations(:,3))
max(correlations(:,3))
histogram(correlations(:,3))

%theta(1) = learning rate
mean(correlations(:,4)) %correlated -.43, max = -.93, min = -.008; left skewed;
min(correlations(:,4))
max(correlations(:,4))
histogram(correlations(:,4))

mean(correlations(:,5))
min(correlations(:,5))
max(correlations(:,5))
histogram(correlations(:,5))

