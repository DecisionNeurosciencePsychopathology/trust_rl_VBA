function  [fx] = f_trust_Qlearn1(x,theta,u,inF)
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
%   and evolution parameter, respectively.F 

r = 1; % when trustee shares the reward is $1.5, or for simplicity r = 1


% theta(1) -- basic learning rate
% theta(2) -- punishment sensitivity

if inF.assymetry_choices==1
%     alpha=1./(1+exp(-theta(1).*u(1)+theta(2).*(u(1)-1)));    
    alpha=1./(1+exp(-theta(1)+theta(2).*(u(1)-1)));    

else
    alpha = 1./(1+exp(-theta(1))); % learning rate is bounded between 0 and 1.
end
% fx = 0;
pe = u(2).*r-x(1); % prediction error

%% introduce reputation sensitivity
if inF.reputation_sensitive==1
    theta(2) = sig(theta(2));
    fx(1) = x(1)+alpha*pe + theta(2).*u(3).*u(4);
elseif inF.humanity==1
    theta(2) = sig(theta(2));
    fx(1) = x(1)+alpha*pe + theta(2).*u(5).*u(4);
elseif inF.valence_n==1 && inF.valence_p==1
%    fx(1) = x(1)+alpha*pe +theta(2).*u(6).*u(4)+theta(3)*u(7).*u(4);
    theta(2) = sig(theta(2));
    theta(3) = sig(theta(3));
    fx(1) = x(1)+alpha*pe -theta(2).*u(7).*u(4)+theta(3).*u(6).*u(4);
elseif inF.valence_n==1
%    fx(1) = x(1)+alpha*pe +theta(2).*u(6).*u(4);
    theta(2) = sig(theta(2));
    fx(1) = x(1)+alpha*pe -theta(2).*u(7).*u(4);
elseif inF.valence_p==1
%    fx(1) = x(1)+alpha*pe +theta(2).*u(7).*u(4);
    theta(2) = sig(theta(2));
    fx(1) = x(1)+alpha*pe +theta(2).*u(6).*u(4);   
else
    fx(1) = x(1) + alpha*pe;
end
%tracking PEs
fx(2) = pe;

% fx(2) = pe;
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


% dfdx = [1-alpha];
% %             0, 0];
% dfdP = [alpha*(1-alpha)*pe];
