function  [ gx] = g_trust_softmax_ED(x,phi,u,in )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(phi(1));
kappa = phi(2);


dQ = x(1);
%dQ = x(1)-1; %for the corrected null model


% dQ = (x(1)-x(2));
gx = sig(kappa + beta*dQ);
% dgdx = zeros(size(x,1),1); % for value tracking only
% % dgdx = zeros(size(1,1),1);  % for value + pe tracking
% dgdx(1) = beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
