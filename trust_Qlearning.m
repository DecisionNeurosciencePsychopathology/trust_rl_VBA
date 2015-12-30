% Q-learning demo
% In psychological terms, motivation can be defined as the set of processes
% that generate goals and thus determine behaviour. A goal is nothing else
% than a ???state of affairs???, to which people attribute (subjective) value.
% Empirically speaking, one can access these values by many means,
% including subjective verbal report or decision making. These measures
% have been used to demonstrate how value change as people learn a new
% operant response. This is predicted by reinforcement learning theories,
% which essentially relate behavioural response frequency to reward. In
% this context, value is expected reward, and it changes in proportion to
% the agent's prediction error, i.e. the difference between actual and
% expected reward.
% This demo simulates a number of sequences of choices of a Q-learning
% agent, which is a simple example of reinforcement learning algiorithms.
% We then invert the model using VBA. Finally, we perform a Volterra
% decomposition of hidden states dynamics onto a set of appropriately
% chosen basis functions (here: the agent's chosen action, and the winning
% action). This diagnostic analysis allows one to identify the hidden
% states' impulse response to experimentally controlled inputs to the
% system.

close all
clear variables




f_fname = @f_trust_Qlearn1; % evolution function (Q-learning)
g_fname = @g_trust_softmax1; % observation function (softmax mapping)
%h_fname = @h_randOutcome; % feedback function (reward schedule)
h_fname = @h_Id; % feedback function, reads from u

multisession = 1;
fixed = 1;
ntrials = 192;

n_hidden_states = 1; %only track the value of sharing, i.e. V(trustee)
%% load subject's data
cd /Users/localadmin/Google' Drive'/skinner/trust/scan_behavior/
%cd /Users/localadmin/trust_rl/data
% our favorite Qlearning subject: 881105
load trust881075.mat
actions = share(1:ntrials)'; %subject's actions
u(1,:) = double(actions);
rewards = double(strcmp(b.TrusteeDecides(b.Order_RS>-999),'share')); %rewards including counterfactual ones (trustee's actions)
rewards(rewards==0) = -1;
u(2,:) = rewards(1:ntrials)';
in.u = u;
y = u(1,:); %the subject's actions

%% allocate feedback struture for simulations
% fb.inH.er = 1;
% fb.inH.vr = 0;
% fb.h_fname = h_fname;
% fb.indy = 1;
% fb.indfb = 2;
% u0 = [randn(1,25)>-0.25]; % possible feedbacks
% fb.inH.u0 = u(2,:); % reinforcement history

%% dimensions and options
% simulation parameters
% theta = sigm(0.75,struct('INV',1)); % learning rate = 0.75
% phi = log(2); % inverse temperature = 2
x0 = zeros(n_hidden_states,1);
% n_t = size(fb.inH.u0,2)+1; % number of trials
n_t = size(u(2,:),2); % number of trials
n_trustees = 4;
if multisession
    options.multisession.split = repmat(n_t/n_trustees,1,n_trustees); % two sessions of 120 datapoints each
    %% fix parameters
    if fixed
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        
        % allow unique initial values for each trustee?
        options.multisession.fixed.X0 = 'all';
    end
    
    
end
%  options.isYout(u(2,:)==-999)=1;
options.binomial = 1;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.
% [y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,zeros(2,n_t),Inf,Inf,options,x0,fb);
% hf = figure('color',[1 1 1]);
% ha = axes('parent',hf,'nextplot','add');
% plot(ha,y,'kx')
% plot(ha,y-e,'r')
% legend(ha,{'y: agent''s choices','p(y=1|theta,phi,m): behavioural tendency'})
%



%% defined number of hidden states and parameters
dim = struct('n',n_hidden_states,'n_theta',1,'n_phi',1, 'n_t', n_t);


%% priors
priors.muPhi = ones(dim.n_phi,1);
priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(n_hidden_states,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
priors.SigmaX0 = .5*eye(dim.n);
priors.a_sigma = 1;       % Jeffrey's prior
priors.b_sigma = 1;       % Jeffrey's prior
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;

options.verbose=1;
options.DisplayWin=1;

%% model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
% displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);
%% select fixed-effects fitting of multi-session data

%
% if ~multisession
% % perform Volterra decomposition
% u1 = u(1,:); % own action
% u3 = u(2,:); % feedback
% u2 = zeros(size(u1)); % opponent's action
% u2(u3>0) = u1(u3>0);
% u2(u3<0) = 1-u1(u3<0);
% uu = 2*[u1;u2]-1;
% o = out;
% o.u = uu;
% [kernels] = VBA_VolterraKernels(posterior,o,16);
% o.diagnostics.kernels = kernels;
% VBA_ReDisplay(posterior,o,1);
%
% getSubplots
% end