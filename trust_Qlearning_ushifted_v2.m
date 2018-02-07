function [posterior,out] = trust_Qlearning_ushifted_v2(datalocation,id, counter, multisession, fixed, sigmakappa, censor, save_str)

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
    counter = 0;
    multisession = 0;
    fixed = 1;    
    sigmakappa = 0;
    censor = 0;
    save_str = '';
elseif nargin<4
    multisession = 0;
    fixed = 1;
    sigmakappa = 0;
    censor = 0;
    save_str = '';
elseif nargin<5
    fixed = 1;
    sigmakappa = 0;
    censor = 0;
    save_str = '';
elseif nargin<6
    sigmakappa = 0;
    censor = 0;
    save_str = '';
elseif nargin<7
    censor = 0;
    save_str = '';
elseif nargin<8
    save_str = '';
end

close all
% clear variables



%% Evolution and observation functions
if counter == 0
    %f_fname = @f_trust_Qlearn1;                        % evolution function (Q-learning) with a single hidden state, Q(share)
    f_fname = @f_trust_Qlearn_null_pmv;                 % new null function that only updates the value of share when the subject shared and has a different choice rule (please note that this latter part needs to be automatic);
elseif counter == 1
    %f_fname = @f_trust_Qlearn_counter_hybrid_regret;   %regret
    f_fname = @f_trust_Qlearn_counter_hybrid_regret_pmv;%regret corrected by PMV    
elseif counter == 2
   %f_fname = @f_trust_Qlearn_mixed_null_countersubject;% evolution function (Q-learning) with a single hidden state, Q(share) and mixed counterfactual + actual rewards
    f_fname = @f_trust_Qlearn_policy;                   % policy: evolution function (Q-learning) with a single hidden state, Q(share) and counterfactual rewards
elseif counter == 3
    f_fname = @f_trust_Qlearn_counter_trustee;          % counterfactual for trustees  
elseif counter == 4
    f_fname = @f_trust_Qlearn_counter_regret2;          % modeling regret based on reviewers' suggestions
else 
    f_fname = @f_trust_Qlearn_counter_disappoint;
 
end

if sigmakappa == 0
    g_fname=@g_trust_softmax1;      % observation function (softmax mapping), evaluates Q(share)
elseif sigmakappa == 1 && counter == 0
    g_fname = @g_trust_softmax_ED1;  % observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for actual evolution function
elseif sigmakappa == 1
     g_fname = @g_trust_softmax_ED2;  % observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter, w/ choice rule appropriate for other than actual evolution function(s) 
else %sigmakappa == 2
    g_fname = @g_trust_softmax_ks_kt;  % observation function (softmax mapping), evaluates Q(share), w/ subject-wise kappa parameter + trustee-wise kappa parameter (requires multisession)
end


%h_fname = @h_randOutcome; % feedback function (reward schedule)
% h_fname = @h_Id; % feedback function, reads from u

n_hidden_states = 2; %track the value of sharing and PE
%n_hidden_states = 1; %only track the value of sharing, i.e. V(trustee)

%% Load subject's data
if strcmp('46069', id)
    id = '046069';
end
load([datalocation sprintf('trust%d',id)])

ntrials = length(b.TrialNumber);

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
rewards = double(strcmp(b.TrusteeDecides,'share')); %rewards including counterfactual ones (trustee's actions)
rewards(rewards==0) = -1;
%rewards(noresponse'==1) = -999;
u(2,:) = rewards(1:ntrials)';

%% Experimental design
%reputation vector
trustee_ID = zeros(length(b.identity),1);
trustee_ID(strcmp(b.identity,'good')) = 1;  % Lina's initial coding = 2
trustee_ID(strcmp(b.identity,'bad')) = -1; 
trustee_ID(strcmp(b.identity,'neutral')) = 0; % Lina: -1
u(3,:) =  trustee_ID(1:ntrials);

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

if ntrials > 192 %include practice trials
   u = u(:,ntrials - 191:end);
   noresponse = noresponse(ntrials - 191:end);
   b.identity = b.identity(ntrials - 191:end); 
end

y = u(1,:); %the subject's actions

%% Optional censoring a block of trials
if censor == 1
    u = u(:,~strcmp(b.identity,'computer'));
    y = y(:,~strcmp(b.identity,'computer'));
    noresponse = noresponse(~strcmp(b.identity,'computer'));
end

%% shifting U
u = [zeros(7,1) u(:,1:end-1)];

%% Sensitivities or bias parameters
% options.inF.reputation_sensitive = reputation_sensitive;

if sigmakappa == 1
    options.inF.kappa_s = 1;        %subject-wise kappa only
elseif sigmakappa == 2
    options.inF.kappa_s = 1;        %subject-wise kappa
    options.inF.kappa_st = 1;       %subject-wise trustee-wise kappa
end

%% allocate feedback structure for simulations
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
% x0 = zeros(n_hidden_states,1);
% n_t = size(fb.inH.u0,2)+1; % number of trials

n_t = size(u(2,:),2); % number of trials
n_trustees = length(unique(b.identity)); 

if censor > 0
    n_trustees = length(unique(b.identity(~strcmp(b.identity,'computer'))));
end

if multisession
    options.multisession.split = repmat(n_t/n_trustees,1,n_trustees); %splitting the sessions
    %% fix parameters
    if fixed == 1
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        %options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        %options.multisession.fixed.X0 = 'all';
        %% set priors on initial states
        priors.SigmaX0 = diag([.3 0]);  %X0 is allowed to vary between sessions
    elseif fixed == 2
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions; subject-wise kappa varies between sessions
        %options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        options.multisession.fixed.X0 = 'all';
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
    elseif fixed == 3
        %options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        %options.multisession.fixed.phi = 'all';
        options.multisession.fixed.theta = 'all';
        %options.multisession.fixed.phi = 1;    %fixing the beta parameter to be the same across sessions
        options.multisession.fixed.phi = 1:2;   %fixing the beta and kappa (subject-specific bias) parameters to be the same across all sessions
        options.multisession.fixed.X0 = 'all';     
        priors.SigmaX0 = diag([0 0]);   %infinite precision priors set on the initial value and PEs
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
n_theta = 1; %evolution function paramaters: by default = 1 (learning rate), adds 1 for each additional parametrized experimental design element

if sigmakappa > 0
    n_phi = 1+sigmakappa;  %observation function parameters: by default = 1 (beta), +1 kappa for the subject's altruistic bias, +1 for the subject's trustee specific bias
else
    n_phi = 1;
end

dim = struct('n',n_hidden_states,'n_theta',n_theta,'n_phi',n_phi, 'n_t', n_t);


%% priors
if n_phi == 2
    priors.muPhi = [1;0];
elseif n_phi == 3
    priors.muPhi = [1;0;0];
else 
    priors.muPhi = 1;
end

priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(n_hidden_states,1);
%if reputation_sensitive||humanity||valence_p||valence_n
%    priors.SigmaTheta = 1e1*eye(dim.n_theta);     
%    priors.SigmaTheta = diag([10 1/3*ones(1,dim.n_theta-1)]);
%elseif assymetry_choices 
%    priors.SigmaTheta = diag([10 10]);
% elseif regret
%     priors.SigmaTheta = diag([10 1/3*ones(1,dim.n_theta-1)]);
% elseif counter == 2
%     priors.SigmaTheta = diag([10 10]);
%else
priors.SigmaTheta = 10;
%end

if sigmakappa > 0
    priors.SigmaPhi = diag([10 1/3*ones(1,dim.n_phi-1)]); %better than when xcompared with sigma_kappa = 1
    %priors.SigmaPhi = diag([10 ones(1,dim.n_phi-1)]);
else
    priors.SigmaPhi = 10;
end



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
%VBA_ReDisplay(posterior,out);
%displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);

%% print condition order for interpreting results
ConditionOrder  = unique(b.identity,'stable');
if censor == 1
    ConditionOrder = ConditionOrder(~strcmp(ConditionOrder, 'computer'));
end

out.design = ConditionOrder;

h = figure(1);
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

filename = [save_str subdir filesep sprintf('%d_cntr%d_mltrun%d_kappa%d_censor%d', id, counter,multisession, sigmakappa, censor)];
figure_filename = [save_str subdir filesep sprintf('%d_cntr%d_mltrun%d_kappa%d_censor%d', id, counter,multisession, sigmakappa, censor)];
savefig(h,figure_filename);
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