function  [ gx,dgdx ] = g_trust_softmax_nobeta(x,phi,u,in )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

%beta = exp(phi(1));
%kappa = phi(2);

dQ = x(1);
% dQ = (x(1)-x(2));
% gx = sig(kappa + beta*dQ);
gx = sig(dQ);
%dgdx = zeros(size(x,1),1);
%dgdx(1) = beta*gx*(1-gx);
% dgdx(2) = -beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
