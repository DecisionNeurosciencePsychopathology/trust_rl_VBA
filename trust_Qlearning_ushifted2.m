function [posterior,out] = trust_Qlearning_ushifted2(id, counter, multisession, fixed, sigmakappa,reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret, SVM)

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
    regret = 0;
    SVM = 0;
elseif nargin<4
    fixed = 1;
    reputation_sensitive = 0;
    sigmakappa = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin<5
    sigmakappa = 0;
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    reputation_sensitive = 0;
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin<6
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    reputation_sensitive = 0;
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin<7
    humanity = 0;
    valence_p = 0;
    valence_n = 0;
    assymetry_choices = 0;
elseif nargin<8
    valence_p =0;
    valence_n = 0;
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin<9
    valence_n = 0;
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin<10
    assymetry_choices = 0;
    regret = 0;
    SVM = 0;
elseif nargin < 11
    regret = 0;
    SVM = 0;
elseif nargin < 12
    SVM = 0;
end

close all
% clear variables

SVM = 1;

%% Evolution and observation functions
if counter == 0
    f_fname = @f_trust_Qlearn1; % evolution function (Q-learning) with a single hidden state, Q(share)
else
    f_fname = @f_trust_Qlearn_counter_hybrid;% evolution function (Q-learning) with a single hidden state, Q(share), and counterfactual rewards
end

if assymetry_choices == 1 || sigmakappa == 0
    g_fname=@g_trust_softmax1;      % observation function (softmax mapping), evaluates Q(share)
else
    g_fname = @g_trust_softmax_ED;  % observation function (softmax mapping), evaluates Q(share), w/ kappa parameter
end

if SVM == 1 
    if counter == 0
        f_fname = @f_trust_SVM1;
        g_fname = @g_trust_SVM1;
    end
end

%h_fname = @h_randOutcome; % feedback function (reward schedule)
% h_fname = @h_Id; % feedback function, reads from u

n_hidden_states = 2; %track the value of sharing and PE
%n_hidden_states = 1; %only track the value of sharing, i.e. V(trustee)


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

%% Load subject's data
%cd /Users/localadmin/Google' Drive'/skinner/trust/scan_behavior/
% exemplar Qlearning subjects: 881105, 213704, 216806, 220043
%cd(datalocation);
if strcmp('46069', id)
    id = '046069';
end
load([datalocation sprintf('trust%s',id)])
%ntrials = length(b.Reversal);
% our favorite Qlearning subject: 881105
% load trust881105.mat;
if exist('filename')
        clear 'filename';
    elseif exist('data_dir_str')
        clear 'data_dir_str';
    elseif exist('stringid')
        clear 'stringid';
    elseif exist('subdir')
        clear 'subdir';
end
if not(exist('decisions'))
    share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
    keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
    missed = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
    b.decisions = zeros(ntrials, 1);
    b.decisions(share) = 1;
    b.decisions(keep) = -1;
    b.decisions(missed) = 0;
end
if exist('share')
    actions = share(1:ntrials)'; %subject's actions
else
    share = zeros(length(b.decisions),1);
    share(b.decisions==1) = 1;
    actions = share(1:ntrials);
end
if not(exist('noresponse'))
    noresponse = zeros(length(b.decisions),1);
    noresponse(b.decisions==0|b.decisions==-999)=1;
end
actions = double(actions);
u(1,:) = actions;

%rewards = double(strcmp(b.TrusteeDecides(b.Order_RS>-999),'share')); %rewards including counterfactual ones (trustee's actions)
actions_trustee = double(strcmp(b.TrusteeDecides,'share')); %trustee shared = 1, trustee kept = 0
%actions_trustee(actions_trustee==0) = -1;
%rewards(noresponse'==1) = -999;
u(2,:) = actions_trustee(1:ntrials)';

%% Experimental design
%reputation vector
trustee_ID = zeros(length(b.identity),1);
trustee_ID(strcmp(b.identity,'good')) = 1;  % Lina's initial coding = 2
trustee_ID(strcmp(b.identity,'bad')) = -1; 
trustee_ID(strcmp(b.identity,'neutral')) = 0; % Lina: -1
u(3,:) =  trustee_ID(1:ntrials);

%trustee prerating
trustee_prerate = zeros(length(b.identity),1);
trustee_prerate(strcmp(b.identity,'good')) = b.rating_pre.good;
trustee_prerate(strcmp(b.identity,'bad')) = b.rating_pre.bad;
trustee_prerate(strcmp(b.identity,'neutral')) = b.rating_pre.neutral;
trustee_prerate(strcmp(b.identity,'computer')) = b.rating_pre.computer;
u(8,:)=trustee_prerate;

% humanity vector
human = ones(length(b.identity),1);
human(strcmp(b.identity,'computer')) = 0;
u(5,:) = human(1:ntrials);

% positive/negative valence vector
valence_pos = zeros(length(b.identity),1);
valence_neg = zeros(length(b.identity),1);
valence_pos(strcmp(b.identity,'good')) = 1;
valence_neg(strcmp(b.identity,'bad')) = 1;
u(6,:) = valence_pos(1:ntrials);
u(7,:) = valence_neg(1:ntrials);

%initial state value only to be sensitive
index_vector = zeros(length(b.identity),1);
blocklength = 48;
index_vector([1,1+blocklength,1+2*blocklength, 1+3*blocklength],1) = 1;
u(4,:) = index_vector(1:ntrials);


y = u(1,:); %the subject's actions
%shifting u
u = [zeros(8,1) u(:,1:end-1)];


%% Sensitivities or bias parameters
options.inF.reputation_sensitive = reputation_sensitive;
options.inF.humanity = humanity;
options.inF.valence_p = valence_p;
options.inF.valence_n = valence_n;
options.inF.assymetry_choices = assymetry_choices;
options.inF.regret = regret;
options.inF.svm_par = SVM;

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
%x0 = zeros(n_hidden_states,1);
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

%% defined number of hidden states and parameters
n_theta = 1+reputation_sensitive+humanity+valence_p+valence_n+assymetry_choices+regret+SVM; %evolution function paramaters: by default = 1 (learning rate), adds 1 for each additional experimental design element
n_phi = 1+sigmakappa; %observation function parameters: by default = 1 (beta), adds 1 for each additional observation parameter
dim = struct('n',n_hidden_states,'n_theta',n_theta,'n_phi',n_phi, 'n_t', n_t);


%% priors
if n_phi ==2
    priors.muPhi = [1;0]; 
else 
    priors.muPhi = 0;           %Fareri et al., 2015: rate of exploration (for choice rule) should vary between 0 and 1
end

%priors.muTheta = zeros(dim.n_theta,1);
%priors.muX0 = zeros(n_hidden_states,1);

priors.muTheta = [0; 0];                %Fareri et al., 2015: learning parameter varies from 0 to 1; SVM_par varies from 0 to 5
priors.muX0 = [0.5*1.5; 0.5];           %Fareri et al., 2015: initial value of sharing and initial value of p(trustee_sharing) = 0.5

%% setting variance learning parameter and SVM for the observation function
if SVM
    priors.SigmaTheta = diag([1 1]);%Fareri et al., 2015: learning parameter varies from 0 to 1; SVM_par varies from 0 to 5
end

%% set priors on initial states
priors.SigmaX0 = diag([0 0]);          %Fareri et al., 2015: tracking expected value and probability of trustee sharing
%priors.SigmaX0 = .3*eye(dim.n);% tracking a single hidden state (value)
%priors.SigmaX0 = diag([.3 0]);  % tracking value and prediction error
%priors.SigmaX0 = 0*eye(dim.n); %used this before for tacking value and
%prediction error
% priors.SigmaX0 = zeros(dim.n); % because of correlation with learning
% % rates, one may want to fix initial states to analyze trustee-wise
% % learning rates

%% set hyper-priors
priors.a_sigma = 1;       % Jeffrey's prior
priors.b_sigma = 1;       % Jeffrey's prior
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;

options.verbose=1;
options.DisplayWin=1;
options.GnFigs=1;

%% model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
%displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);

%% print condition order for interpreting results
ConditionOrder  = unique(b.identity,'stable');
out.design = ConditionOrder;

if not(isnumeric(id))
    id = str2double(id);
end

subdir = func2str(f_fname);

if ~exist(subdir,'dir')
    mkdir(subdir)
end

h = figure(1);

%% get prediction errors
alpha = 1./(1+exp(-posterior.muTheta(1)));
if ~multisession
    out.suffStat.PE = diff(out.suffStat.muX)./alpha;
else
    out.suffStat.PE = diff(sum(out.suffStat.muX))./alpha;
end

filename = [subdir filesep sprintf('%d_cntr%d_SVM%d', id, counter, SVM)];
savefig(h,[subdir filesep sprintf('%d_cntr%d_SVM%d', id, counter, SVM)]);
save(filename); 


%% get prediction error time course
 %[fx,dfdx,dfdP,pe] = f_trust_Qlearn_counter(x,P,u,in)

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