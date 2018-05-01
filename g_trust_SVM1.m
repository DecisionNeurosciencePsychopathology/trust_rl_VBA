function  [ gx] = g_trust_SVM1(x,phi,u,in )
% Alex -- consider adding derivatives later for speed/precision

% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = 1/exp(phi(1));     %let beta be inverse temperature so that you could take advantage of the exponential transform

dQ = x(1);              %should be expected value of participant sharing
gx = (exp(dQ*beta))/(exp(dQ*beta) + exp(beta));

% dgdx = zeros(size(x,1),1); % for value tracking only
% dgdx = zeros(size(1,1),1);  % for value + pe tracking
% dgdx(1) = beta*gx*(1-gx);
% dgdP(1) = beta*dQ*gx*(1-gx);
