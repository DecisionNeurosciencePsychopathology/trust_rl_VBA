function [posterior,out] = trust_Qlearning_ushifted(id, counter, multisession, fixed, sigmakappa,reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices)

% VBA fitting of Qlearning to trust data, using VBA-toolbox
%
% [posterior,out] = trust_Qlearning(id, multisession, fixed)
%
% inputs:
% id = 6-digit ID
% multisession = 1 (default)    treat each trustee as a different run, if 0 all runs are concatenated
% fixed = 1 (default)           fix learning rate and inv. temperature, BUT NOT X0, across trustees
% 
% outputs:
% posterior                     posterior distributions
% out                           fit statistics, diagnostics

if nargin<1
    error('*** Enter 6-digit subject ID ***')
elseif nargin<2
    error('*** Specify the f function type: counterfactuals or not ***')
elseif nargin<3
    multisession = 0;
    fixed = 1;    
    reputation_sensitive = 0;
    sigmakappa = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<4
    fixed = 1;
    reputation_sensitive = 0;
    sigmakappa = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<5
    sigmakappa = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    reputation_sensitive = 0;
    assymetry_choices = 0;
elseif nargin<6
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    reputation_sensitive = 0;
    assymetry_choices = 0;
elseif nargin<7
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<8
    valence_p =0;
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<9
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<10
    assymetry_choices = 0;
end

close all
% clear variables



%% Evolution and observation functions
if counter == 0
    f_fname = @f_trust_Qlearn1; % evolution function (Q-learning) with a single hidden state, Q(share)
else
    f_fname = @f_trust_Qlearn_counter;% evolution function (Q-learning) with a single hidden state, Q(share), and counterfactual rewards
end

if assymetry_choices ==1
    g_fname=@g_trust_softmax1;      % observation function (softmax mapping), evaluates Q(share)
else
    g_fname = @g_trust_softmax_ED;  % observation function (softmax mapping), evaluates Q(share), w/ kappa parameter
end


%h_fname = @h_randOutcome; % feedback function (reward schedule)
% h_fname = @h_Id; % feedback function, reads from u

%n_hidden_states = 2; %track the value of sharing and PE
n_hidden_states = 1; %only track the value of sharing, i.e. V(trustee)

%% Where to look for data
%Quick username check, and path setting
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me, 'polinavanyukov') == 1
        datalocation = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/';
    else
        datalocation = '/Users/localadmin/Google Drive/skinner/trust/scan_behavior/';
    end
ntrials = 192;
%ntrials = 190;
% if length(b.TrusteeDecides)<192
%     b.TrusteeDecides = [b.TrusteeDecides; num2cell(ones(ntrials-length(b.TrusteeDecides))*-999)];
% end


%% Load subject's data
%cd /Users/localadmin/Google' Drive'/skinner/trust/scan_behavior/
% exemplar Qlearning subjects: 881105, 213704, 216806, 220043
cd(datalocation);
load(sprintf('trust%s',id))

% our favorite Qlearning subject: 881105
% load trust881105.mat;

actions = share(1:ntrials)'; %subject's actions
actions = double(actions);
%actions(noresponse'==1)= -999;
u(1,:) = actions;

rewards = double(strcmp(b.TrusteeDecides(b.Order_RS>-999),'share')); %rewards including counterfactual ones (trustee's actions)
rewards(rewards==0) = -1;
%rewards(noresponse'==1) = -999;
u(2,:) = rewards(1:ntrials)';

%% Experimental design
%reputation vector
trustee_ID = zeros(length(b.identity),1);
trustee_ID(strcmp(b.identity,'good')) = 1;  % Lina's initial coding = 2
trustee_ID(strcmp(b.identity,'bad')) = -1; 
trustee_ID(strcmp(b.identity,'neutral')) = 0; % Lina: -1
u(3,:) =  trustee_ID;

% humanity vector
human = ones(length(b.identity),1);
human(strcmp(b.identity,'computer')) = 0;
u(5,:) = human;

% positive/negative valence vector
valence_pos = zeros(length(b.identity),1);
valence_neg = zeros(length(b.identity),1);
valence_pos(strcmp(b.identity,'good')) = 1;
valence_neg(strcmp(b.identity,'bad')) = 1;
u(6,:) = valence_pos;
u(7,:) = valence_neg;

%initial state value only to be sensitive
index_vector = zeros(length(b.identity),1);
blocklength = 48;
index_vector([1,1+blocklength,1+2*blocklength, 1+3*blocklength],1) = 1;
u(4,:) = index_vector;

y = u(1,:); %the subject's actions
%u = [zeros(7,1) u(:,2:end)];
%u(2,:) = [0 u(2,1:end-1)];
u = [zeros(7,1) u(:,1:end-1)];
%u(1:7, noresponse'==1) = -999;

%% Sensitivities or bias parameters
options.inF.reputation_sensitive = reputation_sensitive;
options.inF.humanity = humanity;
options.inF.valence_p = valence_p;
options.inF.valence_n = valence_n;
options.inF.assymetry_choices = assymetry_choices;

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
        options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        options.multisession.fixed.phi = 'all';
        options.multisession.fixed.X0 = 'all';
    end
        
end
options.isYout=zeros(size(u(1,:)));
options.isYout(:,1)=1;
options.isYout(noresponse'==1)=1;
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
n_theta = 1+reputation_sensitive+humanity+valence_p+valence_n+assymetry_choices; %by default = 0, adds 1 for each additional experimental design element
n_phi = 1+sigmakappa; %by default = 1, adds 1 for each additional observation parameter
dim = struct('n',n_hidden_states,'n_theta',n_theta,'n_phi',n_phi, 'n_t', n_t);


%% priors
if n_phi ==2
    priors.muPhi = [1 0];
else 
    priors.muPhi = 1;
end

priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(n_hidden_states,1);
if reputation_sensitive||humanity||valence_p||valence_n
%    priors.SigmaTheta = 1e1*eye(dim.n_theta);
    priors.SigmaTheta = diag([10 1/3*ones(1,dim.n_theta-1)]);
elseif assymetry_choices
    priors.SigmaTheta = diag([10 10]);
else
    priors.SigmaTheta = 10;
end
if sigmakappa
    priors.SigmaPhi = diag([10 sigmakappa/3]);
else
     priors.SigmaPhi = 10;
end

%% set priors on initial states
priors.SigmaX0 = .3*eye(dim.n);
%priors.SigmaX0 = 0*eye(dim.n);
% priors.SigmaX0 = zeros(dim.n); % because of correlation with learning
% % rates, one may want to fix initial states to analyze trustee-wise
% % learning rates

%% set hyper-priors
% priors.a_sigma = 1;       % Jeffrey's prior
% priors.b_sigma = 1;       % Jeffrey's prior
% priors.a_alpha = Inf;
% priors.b_alpha = 0;

options.priors = priors;

options.verbose=1;
options.DisplayWin=1;
options.GnFigs=0;

%% model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
%displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);

%% print condition order for interpreting results
ConditionOrder  = unique(b.identity,'stable')
out.design = ConditionOrder;

h = figure(1);
savefig(h,sprintf('%d_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d', id, counter, multisession, fixed, sigmakappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices));

%% get prediction errors
alpha = 1./(1+exp(-posterior.muTheta(1)));
if ~multisession
    out.suffStat.PE = diff(out.suffStat.muX)./alpha;
else
    out.suffStat.PE = diff(sum(out.suffStat.muX))./alpha;
end
 

filename = sprintf('%d_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d', id, counter,multisession, fixed, sigmakappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
save(filename); 

%% get prediction error time course
% [fx,dfdx,dfdP,pe] = f_trust_Qlearn1(x,P,u,in)

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