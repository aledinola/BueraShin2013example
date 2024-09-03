clear
clc
close all
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))
% Buera & Shin (2013) - Financial Frictions and the Persistence of History: A Quantitative Exploration
% The model of Buera & Shin (2013) requires some reworking before
% computing, see the accompanying pdf.
% I use z to denote the markov shock (as that is standard nomenclature of
% VFI Toolkit), which was called e by BS2013.

% BS2013 set
Params.beta=0.904;
% But this results in too much assets such that the eqm interest rate goes
% negative, so I instead just set
Params.beta=0.7; % (I tried 0.85, but still ended up with interest rate of around -0.02, beta=0.8 got interest rate of essentially 0)
% (BS2013 say beta=0.904 was calibrated to target r=0.045; their table 1)
% (I would guess that it is either because they set to max for assets too
% low, or that there is heavy approximation in the value fn at high levels
% of asset; if the original codes were available you would be able to tell)

CreateFigures=1

% A line I needed for running on the Server
%addpath(genpath('./MatlabToolkits/'))


%% Setting
n_d=0;
n_a=1500; % assets
n_ztaupsi=[40,2,2]; % entrepreneurial ability, tax/subsidy distortions, psi (draw new shock or not)
% Note: grid sizes for z, tau and psi are all hardcoded (you cannot change them unless you
% edit the section below that creates the grids and transition probabilities)

%% Parameters

% Preferences
Params.gamma=1.5; % CES utility param
% Params.beta=0.904;

% Production fn
Params.delta=0.06;
Params.alpha=0.33;
Params.upsilon=0.21;

% Entrepreneurial ability shocks
Params.psi=0.894;
Params.eta=4.15;

% Collateral constraint
Params.lambda=1.35;

% Tax/subsidy distortions
Params.tauplus=0.5;
Params.tauminus=-0.15;
Params.q=1.55;

% Initial values for general eqm parameters
% Params.r=0.04;
% Params.w=1;
% I later overwrote these with something closer to what the initial eqm is (just to reduce runtime to find initial eqm)
Params.r=0.02;
Params.w=0.6;


%% Grids and shock process
d_grid=[];
% grid on assets
a_grid=50*linspace(0,1,n_a)'.^3; % the power of 3 puts more points near zero

% The markov shock is psi*(stay where you are)+(1-psi)*(an iid shock)
% Start with the iid, pg 235 of BS2013 explains that
% "The entrepreneurial ability e is assumed to be a truncated and discretized version
% of a Pareto distribution whose probability density if eta*e^-(eta+1) for e>=1. We discretize
% the support of the ability distribution into 40 grid points {e_1,...,e_40}. Denoting the cdf
% of the original Pareto distribution by M(e)=1-e^(-eta), we choose e_1 and e_38 such that M(e_1)=0.633
% and M(e_38)=0.998. Indexing the grid points by j, we construct e_j to be equidistant from j=1,...,38. 
% The two largest values on the grid are given by e_39 and e_40 which satisfy M(e_39)=0.999 and 
% M(e_40)=0.9995. Finally the corresponding probability mass w(e_j) for 2<=j<=40 is given by 
% [M(e_j)-M(e_{j-1})]/M(e_40) and w(e_1)=M(e_1)/M(e_40)."
% Matlab has gpinv to get inverse cdf of generalized pareto distribution, the help give following two lines
% x = gpinv(p,k,sigma,theta): y​​  =​ f(x∣k,σ,θ)=​​​​​​ (1/σ)(1+k(x−θ)/σ)^(−1−1/k)
% When k > 0 and theta = sigma/k, the GP is equivalent to a Pareto distribution with a scale parameter equal to sigma/k and a shape parameter equal to 1/k.
% Pareto dist cdf(x): 1-(scale/x)^shape [got this from wikipedia]
% We want pareto dist cdf as 1-x^(-eta)
% So we want Pareto dist with scale parameter=1, and shape paremeter=eta
% So we want k=1/eta (to hit shape param), sigma=1/eta (to hit scale param), theta=sigma/k=1
zprobs=[linspace(0.633,0.998,38),0.999,0.9995]'; % M(e_j)
z_grid = gpinv(zprobs,1/Params.eta,1/Params.eta,1); % inverse cdf of generalized pareto distribution
% Finally,
pi_z=zeros(40,1);
pi_z(1)=zprobs(1)/zprobs(40);
for jj=2:40
    pi_z(jj)=(zprobs(jj)-zprobs(jj-1))/zprobs(40);
end
% sum(pi_z) % double-check, yes this formula does give total probability of one :)


% Next we have the tax/subsidy distortion shocks
% "we assume that tau is a random variable with only two possible outcomes:
% tauplus (>=0) and tauminus (<=0). Also, the probability of being taxe for
% a type z individual Pr(tau=tauplus|z), is assumed to be 1-exp(-qz).
% Finally, we assume that the idiosyncrati distortions are also governed by
% the same psi shock that determines the persisitence of the
% entrepreneurial productivity, In fact, [individuals] draw a new tau exactly
% when they draw a new z"
tau_grid=[Params.tauplus; Params.tauminus];
% Because the probability of tau depends on the value of z, we cannot create pi_tau, and instead we create pi_taugivenz
pi_taugivenz=zeros(n_ztaupsi(2),n_ztaupsi(1)); % size is tau by z
for z_c=1:n_ztaupsi(1)
    pi_taugivenz(:,z_c)=[1-exp(-Params.q*z_grid(z_c)); exp(-Params.q*z_grid(z_c))]; % Note, second row is just 1-first row 
end

% Put pi_z and pi_taugivenz together to get pi_tauz
pi_ztau=repmat(pi_z,n_ztaupsi(2),1).*reshape(pi_taugivenz',[n_ztaupsi(1)*n_ztaupsi(2),1]);
% (because they are iid (conditional on drawing a new one) I can just do this as a column)

% The third piece is the psi shock that determines whether or not we draw a
% new value. This is needed to impose the 'exactly when' of changing tau
% iff changing z.
psi_grid=[0;1]; % 1 means change z and tau, 0 means fixed
pi_psi=[Params.psi, 1-Params.psi]; % rows are identical as it is iid

% We now need to put together the transition probabilites for (z,tau,psi)
pi_ztaupsi=[pi_psi(1)*eye(n_ztaupsi(1)*n_ztaupsi(2)), pi_psi(2)*pi_ztau'.*ones(n_ztaupsi(1)*n_ztaupsi(2),1); pi_psi(1)*eye(n_ztaupsi(1)*n_ztaupsi(2)), pi_psi(2)*pi_ztau'.*ones(n_ztaupsi(1)*n_ztaupsi(2),1)];
% And stack the grids
ztaupsi_grid=[z_grid; tau_grid; psi_grid];

%% Return fn
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,z,tau,psi,gamma,w,r,lambda,delta,alpha,upsilon) ...
    BueraShin2013_ReturnFn(aprime,a,z,tau,psi,gamma,w,r,lambda,delta,alpha,upsilon);

%% Create some FnsToEvaluate
FnsToEvaluate.A=@(aprime,a,z,tau,psi) a; % assets
FnsToEvaluate.K=@(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon) BueraShin2013_capitaldemand(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon); % Capital used (zero if worker)
FnsToEvaluate.L=@(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon) BueraShin2013_labordemand(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon); % Labor demand (zero if worker)
FnsToEvaluate.entrepreneur=@(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon) BueraShin2013_entrepreneur(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon);% 1 if entrepreneur, 0 if worker
% FnsToEvaluate.z=@(aprime,a,z,tau,psi) z; % entrepreneurial ability

%% Model is set, we can start by just that the basics are running okay
vfoptions=struct();
vfoptions.verbose=1;
vfoptions.lowmemory=0;
vfoptions.parallel=0;
simoptions=struct();

disp('Test a few things before we start')
tic;
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_ztaupsi,d_grid,a_grid,ztaupsi_grid,pi_ztaupsi,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_ztaupsi,pi_ztaupsi,simoptions);
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);
toc

if CreateFigures==1
    % take a look at cdf over asset grid to make sure not hitting top of grid
    figure(1)
    subplot(3,1,1); plot(a_grid,cumsum(sum(sum(sum(StationaryDist,4),3),2)))
    title('cdf of asset to make sure grid on assets seems okay (pre test)')
end
% before we do the general eqm, just take a look at some things to get a
% feel for what going to happen with general eqm conditions
[AggVars.A.Mean,AggVars.K.Mean]
[AggVars.L.Mean,1-AggVars.entrepreneur.Mean]


%% Set up general equilibrium
GEPriceParamNames={'r','w'};

GeneralEqmEqns.capitalmarket=@(K,A) A-K; % assets minus capital demand
GeneralEqmEqns.labormarket=@(L,entrepreneur) L-(1-entrepreneur); % labor demand=labor supply, suppy is just fraction of workers (who each exogneously supply endowment 1 of labor)

%% Now compute the initial stationary general eqm
heteroagentoptions.verbose=1;
vfoptions.verbose=0;
disp('Solving initial stationary general eqm')
tic;
[p_eqm_init,~,GenEqmConds_init]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_ztaupsi, 0, pi_ztaupsi, d_grid, a_grid, ztaupsi_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
toc

% Compute some things relating to initial stationary general eqm
Params.r=p_eqm_init.r;
Params.w=p_eqm_init.w;
[V_init,Policy_init]=ValueFnIter_Case1(n_d,n_a,n_ztaupsi,d_grid,a_grid,ztaupsi_grid,pi_ztaupsi,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_init=StationaryDist_Case1(Policy_init,n_d,n_a,n_ztaupsi,pi_ztaupsi,simoptions);
AggVars_init=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_init, Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy_init, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);

% We will need dist for the transition path
% I want to keep the worker/entrepreneur decision for a graph
workerORentrepreneur_init=ValuesOnGrid.entrepreneur;
clear ValuesOnGrid

if CreateFigures==1
    % take another a look at cdf over asset grid to make sure not hitting top of grid
    figure(1)
    subplot(3,1,2); plot(a_grid,cumsum(sum(sum(sum(StationaryDist_init,4),3),2)))
    title('cdf of asset to make sure grid on assets seems okay (init eqm)')
end


%% And the final stationary general eqm

% Reform eliminates all tau. Easiest way to implement this is just set the
% values of tau_grid to zero (probabilities are then irrelevant; we could
% do 'better' by just eliminating tau from model, but that would be more
% work ;)
tau_grid=[0;0];
ztaupsi_grid=[z_grid; tau_grid; psi_grid];

% % Switch to using shooting algorithm
% heteroagentoptions.fminalgo=5;
% % Need to explain to heteroagentoptions how to use the GeneralEqmEqns to update the general eqm prices.
% heteroagentoptions.fminalgo5.howtoupdate={...  % a row is: GEcondn, price, add, factor
%     'capitalmarket','r',0,0.03;...  % capitalmarket GE condition will be positive if r is too big, so subtract
%     'labormarket','w',1,0.1;... % labormarket GE condition will be positive if w is too small, so add
%     };

disp('Solving final stationary general eqm')
tic;
[p_eqm_final,~,GenEqmConds_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_ztaupsi, 0, pi_ztaupsi, d_grid, a_grid, ztaupsi_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
toc

% Compute some things relating to initial stationary general eqm
Params.r=p_eqm_final.r;
Params.w=p_eqm_final.w;
[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_ztaupsi,d_grid,a_grid,ztaupsi_grid,pi_ztaupsi,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_ztaupsi,pi_ztaupsi,simoptions);
AggVars_final=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_final, Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy_final, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);

% We will need V_final for the transition path
% I want to keep the worker/entrepreneur decision for a graph
workerORentrepreneur_final=ValuesOnGrid.entrepreneur;
clear ValuesOnGrid

if CreateFigures==1
    % take another a look at cdf over asset grid to make sure not hitting top of grid
    figure(1)
    subplot(3,1,3); plot(a_grid,cumsum(sum(sum(sum(StationaryDist_final,4),3),2)))
    title('cdf of asset to make sure grid on assets seems okay (final eqm)')
end

%% This is just me making sure I understand what r and w do to the general eqm conditions
% % Turns out when r goes negative things get messy
% Params.w=0.93;
% for currr=-0.03:0.01:0.03
%     Params.r=currr;
%     [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_ztaupsi,d_grid,a_grid,ztaupsi_grid,pi_ztaupsi,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
%     StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_ztaupsi,pi_ztaupsi,simoptions);
%     AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);
%     Params.r
%     [AggVars.A.Mean,AggVars.K.Mean]
%     [AggVars.L.Mean,1-AggVars.entrepreneur.Mean]
% end
% % r=[-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03]
% % A=[0.9610, 0.9458, 0.9355, 0.9293, 0.9328, 0.9416, 0.9529] % so increasing r increase A, but only for positive r
% % K=[1.2691, 1.2304, 1.1655, 1.0907, 0.9399, 0.9194, 0.0905] % so increasing r decreases K, monotonically
% % When r increases, A increases (as long as r is positive) and K decreases, so A-K will increase
% % [when r inceases, L decreases, 1-entrepreneur increases]
% 
% Params.r=0.02;
% for currw=0.9:0.01:0.96
%     Params.w=currw;
%     [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_ztaupsi,d_grid,a_grid,ztaupsi_grid,pi_ztaupsi,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
%     StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_ztaupsi,pi_ztaupsi,simoptions);
%     AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_ztaupsi, d_grid, a_grid, ztaupsi_grid, [], simoptions);
%     Params.w
%     [AggVars.A.Mean,AggVars.K.Mean]
%     [AggVars.L.Mean,1-AggVars.entrepreneur.Mean]
% end
% %     w=[0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96]
% %     L=[0.9703, 0.9069, 0.8792, 0.8526, 0.8273, 0.8032, 0.7793] % so increasing w decreases L
% % 1-ent=[0.9443, 0.9541, 0.9547, 0.9553, 0.9558, 0.9563, 0.9568] % so increasing w, increases 1-entrepreneur
% % When w increase, L decreases, 1-ent increases, so L-(1-ent) will decrease



%% Before we do the transition path, just a few things that help understand what is going on

% First, a plot of the probability of being taxed (vs subsidised)---the
% 'distortionary wedge'--- as a function of the entrepreneurial abilities
if CreateFigures==1
    figure(2)
    plot(z_grid,pi_taugivenz(1,:),'x')
    title('Probability of tax (vs subsidy) as function of z grid (entrepreneurial ability)')
end
% Shows that high ability entrepreneurs are basically guaranteed to be taxed

% Second, plots before and post-reform (initial and final stationary general eqms)
% of who becomes a worker vs entpreneur (in terms of assets and entrepreneurial ability; a and z)
% Note: You will need to make this figure full screen to be able to read it
% If we plot all the asset grid the ylabels is a mess because of the lines (between each square), so just plot every 10th asset grid point
% (better would plot all and then modify the y-axis so that doesn't ytick them all, but I can't be bothered)
if CreateFigures==1
    just10thasset=1:10:n_a;
    figure(3)
    temp1=gather(workerORentrepreneur_init(just10thasset,:,1,1)); % heatmap only works with cpu
    subplot(3,1,1); heatmap(z_grid,a_grid(just10thasset),temp1)
    grid off
    title('Initial eqm, with tax: Who becomes entrepreneur')
    temp1=gather(workerORentrepreneur_init(just10thasset,:,2,1)); % heatmap only works with cpu
    subplot(3,1,2); heatmap(z_grid,a_grid(just10thasset),temp1)
    grid off
    title('Initial eqm, without tax: Who becomes entrepreneur')
    temp1=gather(workerORentrepreneur_final(just10thasset,:,1,1));% heatmap only works with cpu
    subplot(3,1,3); heatmap(z_grid,a_grid(just10thasset),temp1)
    grid off
    title('Final eqm: Who becomes entrepreneur')
end
% Can see that for initial eqm with tax, only very high entrepreneurial ability become entrepreneurs, and only those with some assets
% Whereas with subsidy, all entrepreneurial abilities become entrepreneurs, as long as they have some assets
% In final eqm, is more like the top two-thirds of entrepreneurial ability, but with a more important interaction with assets

save ./SavedOutput/BS2013_preTpath.mat


%% And now the transition path
% Ths is a little odd for the toolkit, as really the path is on the exogenous tau shock
% But the toolkit will just take these from the inputs, so we can just put the
% final tau shocks (of zero) into Params, and put an irrelevant placebo into 
% ParamPath (which here is lambda, but could be any parameter that is
% unchanged over the transition)

T=100; % number of periods for transition

% Placebo path (lambda is unchanged between initial and final)
ParamPath.lambda=Params.lambda*ones(1,T);

% Some naive guesses for price path
PricePath0.r=[linspace(p_eqm_init.r,p_eqm_final.r,floor(T/2)),p_eqm_final.r*ones(1,T-floor(T/2))];
PricePath0.w=[linspace(p_eqm_init.w,p_eqm_final.w,floor(T/2)),p_eqm_final.w*ones(1,T-floor(T/2))];

% General eqm conditions are actually just same
TransPathGeneralEqmEqns.capitalmarket=@(K,A) A-K; % capital demand equals assets
TransPathGeneralEqmEqns.labormarket=@(L,entrepreneur) L-(1-entrepreneur); % labor demand=labor supply, suppy is just fraction of workers (who each exogneously supply endowment 1 of labor)

% But now we set them up as shooting-algorithm
transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to
% update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate={... % a row is: GEcondn, price, add, factor
    'capitalmarket','r',0,0.03;...  % capitalmarket GE condition will be positive if r is too big, so subtract
    'labormarket','w',1,0.1;... % labormarket GE condition will be positivie if w is too small, so add
    };

% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.


% Now just run the TransitionPath_Case1 command (all of the other inputs are things we 
% had already had to define to be able to solve for the initial and final equilibria)
transpathoptions.verbose=1;
if CreateFigures==1
    transpathoptions.graphpricepath=1;
    transpathoptions.graphaggvarspath=1;
end

disp('Solving transition path')
tic;
PricePath=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_ztaupsi, pi_ztaupsi, d_grid,a_grid,ztaupsi_grid, ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions);
toc

save ./SavedOutput/BS2013.mat

% Calculate some things about the transition
[VPath,PolicyPath]=ValueFnOnTransPath_Case1(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_ztaupsi, pi_ztaupsi, d_grid, a_grid,ztaupsi_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions, simoptions);
AgentDistPath=AgentDistOnTransPath_Case1(StationaryDist_init, PolicyPath,n_d,n_a,n_ztaupsi,pi_ztaupsi,T,simoptions);

% Figure 3 shows GDP, TFP and Investment rate
% Figure 4 shows Capital Stock and Interest Rates
% Figure 5 shows Average Entrepreneurial Ability and Wealth Share of Top 5% Ability

% Create some additional FnsToEvaluate so we can get all of these
FnsToEvaluate.abilityofentrepreneur=@(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon) z*BueraShin2013_entrepreneur(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon);% 1 if entrepreneur, 0 if worker
% Note: To get the average we will need to renormalize by fraction that are entrepreneurs


AggVarsPath=EvalFnOnTransPath_AggVars_Case1(FnsToEvaluate,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_ztaupsi, pi_ztaupsi, d_grid, a_grid,ztaupsi_grid,simoptions);





















