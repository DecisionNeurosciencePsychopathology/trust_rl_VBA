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


