function  [fx] = f_trust_Qlearn_counter_hybrid(x,theta,u,inF)

% function  [fx,dfdx,dfdP] = f_trust_Qlearn_counter(x,theta,u,inF)
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
%   - dfdx/dfdP: gradient of the q-values evolution function, wrt q-values
%   and evolution parameter, respectively.F 

%r = 1; % when trustee shares the reward is $1.5, or for simplicity r = 1
% theta(1) -- basic learning rate
% theta(2) -- punishment sensitivity or reality of rewards parameter

%write out actual rewards first, then counterfactual w/ omega

% %% code actual rewards
% actual_reward = 0;
% if (u(2)==1 && u(1)==1)     %trustee shared, subject shared
%     actual_reward = 1.5;
% elseif (u(2)==0 && u(1)==1) %trustee kept, subject shared
%     actual_reward = 0;
% elseif (u(2)==0 && u(1)<1)  %trustee kept, subject kept
%     actual_reward = -1;
% else 
%     actual_reward = -1;
% end
% 
% %% counterfactual rewards
% counter_reward = 0;
% if (u(2)==1 && u(1)==1)     %trustee shared, subject shared
%     counter_reward = 0;
% elseif (u(2)==0 && u(1)==1) %trustee kept, subject shared
%     counter_reward = -1;
% elseif (u(2)==0 && u(1)<1)  %trustee kept, subject kept
%     counter_reward = -1;
% else
%     counter_reward = 0.5;
% end

if (u(2)==1 && u(1)==1)     %trustee shared, subject shared
    r = 1.5;                % BUT why don't they consider their alternative action as a reference?
elseif (u(2)<1 && u(1)==1) %trustee kept, subject shared
    r = -1.5;
elseif (u(2)<1 && u(1)<1)  %trustee kept, subject kept
    r = -1;
else r = 0.5;               %trustee shared, subject kept
end

%% code counterfactual rewards (also apply cr to all trials OR incongruent trials)?

% a parameter the modulates the "regret" wrt "share" action
if inF.regret == 1
    %omega = 1./(1+exp(-theta(2))); %bounded between 0 and 1.
    omega1 = theta(2);
%    omega2 = theta(3);
    if (u(2)<1 && u(1) ==1)         %trustee kept, subject shared
        r = actual_reward + counter_reward * (1+omega1);
    elseif  (u(2)==1 && u(1) <1)    %trustee shared, subject kept
        r = actual_reward + counter_reward * (1+omega1);
    else
        r = actual_reward + counter_reward;
    end
else
    %r = actual_reward + counter_reward;
    r = r;
end

if inF.assymetry_choices==1
%     alpha=1./(1+exp(-theta(1).*u(1)+theta(2).*(u(1)-1)));    
    alpha=1./(1+exp(-theta(1)+theta(2).*(u(1)-1)));    
else
    alpha = 1./(1+exp(-theta(1))); % learning rate is bounded between 0 and 1.
end


pe = r-x(1); % prediction error
fx = zeros(length(x),1);

%% introduce reputation sensitivity: this assumes that reputation sensitivity is
%% an additive effect wrt the initial value state

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
    fx(1) = x(1)+alpha*pe + theta(2).*u(6).*u(4);   
else
    fx(1) = x(1) + alpha*pe;
end

%tracking PEs
fx(2) = pe;

%% one hidden state (value)
% dfdx = zeros(size(x,1),1);
% dfdx(1) = [1-alpha];
% dfdP = [alpha*(1-alpha)*pe];

%% two hidden states (value + pe)
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
