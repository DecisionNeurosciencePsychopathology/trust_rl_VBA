
[phi1,phi2]=meshgrid(-2:0.1:2, -2:0.1:2);

dQ = [-.5 .5];
models = {'old' 'new'};
beta = exp(phi1);
kappa = phi2;

is_new_parameterization = 0;

figure(100); clf;
i = 1;

for q = 1:length(dQ)
    dq = dQ(q);
    for model = 1:length(models);
        p_share = learning_rule(phi1,kappa,dq,models(model));
        subplot(2,2,i);
        i = i+1;
        surface(beta, kappa, p_share);
        view(3);
        title(sprintf('%s model, value=%d', char(models(model)), dq));
        fprintf('plotted one  ');
    end
end
