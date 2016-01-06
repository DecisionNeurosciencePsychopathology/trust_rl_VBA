function  [fx,dfdx,dfdP] = f_trust_Qlearn1(x,theta,u,in)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : q-values (2x1)
%   - P : (inverse-sigmoid) learning-rate
%   - u : u(1)=previous action (1 or 0), u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved q-values (2x1)
%   - dfdx/dfdP: gradient of the q-values evolution function, wrt q-avlues
%   and evolution parameter, respectively.
r = 1; % when trustee shares the reward is $1.5, or for simplicity r = 1
alpha = 1./(1+exp(-theta)); % learning rate is bounded between 0 and 1.
% fx = 0;
pe = u(2).*r-x; % prediction error

% update the value of SHARE
fx = x + alpha*pe;
% the reward for KEEP is always $1
% fx(2) = 1;

% gradients' derivation
% if u(1)==1
%     dfdx = [1-alpha, 0;
%             0, 1];
%     dfdP = [alpha*(1-alpha)*pe(1),0];
% else
%     dfdx = [1, 0;
%             0, 1-alpha];
%     dfdP = [0,alpha*(1-alpha)*pe(2)];
% end


dfdx = [1-alpha];
%             0, 0];
dfdP = [alpha*(1-alpha)*pe];
