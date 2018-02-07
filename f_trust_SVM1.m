function  [fx] = f_trust_SVM1(x,theta,u,inF)
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

% theta(1) -- basic learning rate
% theta(2) -- SVM_parameter

if inF.assymetry_choices==1
%     alpha=1./(1+exp(-theta(1).*u(1)+theta(2).*(u(1)-1)));    
    alpha=1./(1+exp(-theta(1)+theta(2).*(u(1)-1)));    

else
    alpha = 1./(1+exp(-theta(1)));   %Fareri et al., 2015: default learning rate is bounded between 0 and 1.
end

    SVM_par = 5./(1+exp(-theta(1))); %Fareri et al., 2015: SV parameter is bounded between 0 and 1 and scaled up by a factor of 5.

pe = u(2)-x(2);                         %where x(2) should be probability of trusteee sharing calculated from previous trials
fx(2) = x(2)  + alpha * pe;             %Fareri et al., 2015: calculating the probability of trustee sharing;
fx(1) = fx(2)*(1.5 + (SVM_par.*u(8)/7));%Fareri et al., 2015: expected value of sharing based on the probability of trustee sharing and SVM parameter;

