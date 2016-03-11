function  [ gx,dgdx, dgdP ] = g_trust_softmax2(x,phi,u,in )
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
% dQ = (x(1)-x(2));
%gx = sig(kappa + beta*dQ);
gx = sig(beta*(kappa+dQ));
dgdx = zeros(size(x,1),1);
dgdx(1) = (beta*(exp(-beta*dQ-beta*kappa)))/((exp(-beta*dQ-beta*kappa)+1)^2); %derivative wrt dQ
dgdP(1) = ((kappa+dQ)*exp(-(kappa+dQ)*beta))/(exp(-(kappa+dQ)*beta) + 1)^2; %derivative wrt dP=beta
dgdP(2) = (beta*(exp(-beta*dQ-beta*kappa)))/((exp(-beta*dQ-beta*kappa)+1)^2); %derivative wrt dP=kappa
% dgdx(2) = -beta*gx*(1-gx);
 %dgdP(1) = beta*gx*(1-gx);
