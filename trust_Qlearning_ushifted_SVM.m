function [posterior,out] = trust_Qlearning_ushifted_SVM(datalocation, id, counter, multisession, fixed, sigmakappa, save_str, SVM)

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
f_fname = @f_trust_SVM1;
g_fname = @g_trust_SVM1;


n_hidden_states = 2; %track the value of sharing and PE
%n_hidden_states = 1; %only track the value of sharing, i.e. V(trustee)

%% Load subject's data
%cd /Users/localadmin/Google' Drive'/skinner/trust/scan_behavior/
% exemplar Qlearning subjects: 881105, 213704, 216806, 220043
%cd(datalocation);
if strcmp('46069', id)
    id = '046069';
end
load([datalocation sprintf('trust%d',id)])

ntrials = length(b.TrialNumber);
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
%trustee_prerate(strcmp(b.identity,'computer')) = b.rating_pre.computer;

%pre-ratings are normalized in the evolution function
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
n_theta = 1+SVM; %evolution function paramaters: learning rate in Eq. 3, scaling or SVM parameter in Eq. 5 (Fareri et al., 2015)
n_phi = 1; %observation function parameters: by default = 1 (beta), adds 1 for each additional observation parameter
dim = struct('n',n_hidden_states,'n_theta',n_theta,'n_phi',n_phi, 'n_t', n_t);


%% priors
priors.muPhi = 1;               %Fareri et al., 2015: rate of exploration (for choice rule) should vary between 0 and 1
priors.muTheta = [0; 0];        %Fareri et al., 2015: learning parameter varies from 0 to 1; SVM_par varies from 0 to 5 (re-scaling occurs inside the evolution function)
priors.muX0 = [0.5*1.5; 0.5];   %Fareri et al., 2015: initial value of sharing and initial value of p(trustee_sharing) = 0.5

priors.SigmaPhi = 10;
priors.SigmaTheta = diag([10 10]);%Fareri et al., 2015: learning parameter varies from 0 to 1; SVM_par varies from 0 to 5
priors.SigmaX0 = diag([0 0]);   %Fareri et al., 2015: tracking expected value and probability of trustee sharing - fixed.


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

%% get prediction errors
alpha = 1./(1+exp(-posterior.muTheta(1)));
if ~multisession
    out.suffStat.PE = diff(out.suffStat.muX)./alpha;
else
    out.suffStat.PE = diff(sum(out.suffStat.muX))./alpha;
end


if not(isnumeric(id))
    id = str2double(id);
end

subdir = func2str(f_fname);

if ~exist([save_str subdir],'dir')
    mkdir([save_str subdir])
end

filename = [save_str subdir filesep sprintf('%d_cntr%d_mltrun%d_fixed%d_kappa%d', id, counter,multisession, fixed, sigmakappa)];
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