

function p_share = learning_rule(phi1,kappa,dQ,model)
%  learning rule for plotting diagnostics
beta=exp(phi1);
if strcmpi(model,'new')
    p_share = sig(beta.*(kappa + dQ));
    
elseif strcmpi(model,'old')
    p_share = sig(kappa + beta.*dQ);
end
