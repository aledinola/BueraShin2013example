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
%Params.beta=0.904;
% But this results in too much assets such that the eqm interest rate goes
% negative, so I instead just set
%Params.beta=0.7; % (I tried 0.85, but still ended up with interest rate of around -0.02, beta=0.8 got interest rate of essentially 0)
% (BS2013 say beta=0.904 was calibrated to target r=0.045; their table 1)
% (I would guess that it is either because they set to max for assets too
% low, or that there is heavy approximation in the value fn at high levels
% of asset; if the original codes were available you would be able to tell)

CreateFigures=0;

%% Setting
n_d = 0;
n_a = 1000; % assets
n_z = 40; % entrepreneurial ability

%% Parameters

% Preferences
Params.crra = 1.5;   % CRRA utility param
Params.beta = 0.904; % Discount factor

% Production fn
Params.delta   =0.06; % Capital depreciation rate
Params.alpha   =0.33; % Capital coefficient
Params.upsilon =0.21; % 1 minus span of control parameter

% Entrepreneurial ability shocks
Params.psi     =0.894; %stochastic process: persistence
Params.eta     =4.15;  %stochastic process: dispersion

% Collateral constraint
Params.lambda  =inf;%1.35;

% Initial values for general eqm parameters
% Params.r=0.04;
% Params.w=1;
% I later overwrote these with something closer to what the initial eqm is (just to reduce runtime to find initial eqm)
Params.r = 0.0472;
Params.w = 0.171;

%% Grid for assets
d_grid=[];
% grid on assets
n_a    = 1001;
a_min  = 1e-6;
a_max  = 4000;
% a_scale>1 puts more points near zero
a_scale = 2;
a_grid  = a_min+(a_max-a_min)*linspace(0,1,n_a)'.^a_scale;

%% Stochastic process
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
% zprobs=[linspace(0.633,0.998,n_z-2),0.999,0.9995]'; % M(e_j)
% z_grid = gpinv(zprobs,1/Params.eta,1/Params.eta,1); % inverse cdf of generalized pareto distribution
% % Finally,
% pi_z_vec=zeros(n_z,1);
% pi_z_vec(1)=zprobs(1)/zprobs(n_z);
% for jj=2:n_z
%     pi_z_vec(jj)=(zprobs(jj)-zprobs(jj-1))/zprobs(n_z);
% end
% % sum(pi_z_vec) % double-check, yes this formula does give total probability of one :)
% if abs(1-sum(pi_z_vec))>1e-8
%     disp(sum(pi_z_vec))
%     error('pi_z_vec does NOT sum to one!')
% end

z_grid   = importdata('support.dat');
pi_z_vec = importdata('dist.dat');
% See my notes on Buera and Shin paper
pi_z = Params.psi*eye(n_z)+(1-Params.psi)*ones(n_z,1)*pi_z_vec';
pi_z = pi_z./sum(pi_z,2);

%% Return fn
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,z,crra,w,r,lambda,delta,alpha,upsilon) ...
    BueraShin2013_ReturnFn(aprime,a,z,crra,w,r,lambda,delta,alpha,upsilon);

%% Create some FnsToEvaluate
FnsToEvaluate.A=@(aprime,a,z) a; % assets
% Capital used (zero if worker)
FnsToEvaluate.K=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_capitaldemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% Labor demand (zero if worker)
FnsToEvaluate.L=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_labordemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% 1 if entrepreneur, 0 if worker
FnsToEvaluate.entrepreneur=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_entrepreneur(aprime,a,z,w,r,lambda,delta,alpha,upsilon);
% Entrepreneurial output (zero if worker)
FnsToEvaluate.Y=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_output(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% Enternal finance
FnsToEvaluate.extfin=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_extfin(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 

%% Model is set, we can start by just that the basics are running okay
vfoptions           = struct();
vfoptions.verbose   = 1;
vfoptions.lowmemory = 0;
simoptions          = struct();

% disp('Test a few things before we start')
% tic;
% [V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
% StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
% AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
% toc
% 
% ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy,FnsToEvaluate,...
%     Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,[],simoptions);
% 
% % 1=entre, 0 = worker
% workerORentrepreneur = gather(ValuesOnGrid.entrepreneur);
% clear ValuesOnGrid
% 
% if CreateFigures==1
%     % take a look at cdf over asset grid to make sure not hitting top of grid
%     figure
%     plot(a_grid,cumsum(sum(StationaryDist,2)))
%     title('cdf of asset to make sure grid on assets seems okay (pre test)')
% 
%     % Occupational choice 1
%     just10thasset=1:10:n_a;
%     figure
%     temp1 = workerORentrepreneur(just10thasset,:); 
%     heatmap(z_grid,a_grid(just10thasset),temp1)
%     grid off
%     title('Who becomes entrepreneur')
%     xlabel('Ability')
%     ylabel('Assets')
% 
% end
% 
% % Before we do the general eqm, just take a look at some things to get a
% % feel for what going to happen with general eqm conditions
% disp('Look at GE conditions: capital, labor')
% disp([AggVars.A.Mean,AggVars.K.Mean])
% disp([AggVars.L.Mean,1-AggVars.entrepreneur.Mean])

%% Set up general equilibrium
GEPriceParamNames={'r','w'};

GeneralEqmEqns.capitalmarket=@(K,A) A-K; % assets minus capital demand
GeneralEqmEqns.labormarket=@(L,entrepreneur) L-(1-entrepreneur); % labor demand=labor supply, suppy is just fraction of workers (who each exogneously supply endowment 1 of labor)

%% Now compute the initial stationary general eqm
heteroagentoptions.verbose=1;
heteroagentoptions.toleranceGEprices=10^(-3);
heteroagentoptions.toleranceGEcondns=10^(-3);
vfoptions.verbose=1; 
do_GE = 0;

if do_GE==1
    disp('Solving initial stationary general eqm')
else
    disp('Solving initial partial eqm')
end

[Outputs,GE_cond,Policy_init,StationaryDist_init,AggVars_init,ValuesOnGrid] = BueraShin_Fn(do_GE,Params,n_d,n_a,n_z,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);

%% Initial stationary general eqm
Params.r=Outputs.r;
Params.w=Outputs.w;
% We will need dist for the transition path
% I want to keep the worker/entrepreneur decision for a graph
% 1=entre, 0=worker
workerORentrepreneur_init=ValuesOnGrid.entrepreneur;
clear ValuesOnGrid

if CreateFigures==1
    % take another a look at cdf over asset grid to make sure not hitting top of grid
    figure(1)
    subplot(3,1,2); plot(a_grid,cumsum(sum(sum(sum(StationaryDist_init,4),3),2)))
    title('cdf of asset to make sure grid on assets seems okay (init eqm)')
end

% % Switch to using shooting algorithm
% heteroagentoptions.fminalgo=5;
% % Need to explain to heteroagentoptions how to use the GeneralEqmEqns to update the general eqm prices.
% heteroagentoptions.fminalgo5.howtoupdate={...  % a row is: GEcondn, price, add, factor
%     'capitalmarket','r',0,0.03;...  % capitalmarket GE condition will be positive if r is too big, so subtract
%     'labormarket','w',1,0.1;... % labormarket GE condition will be positive if w is too small, so add
%     };


%% Analize results

% Second, plots before and post-reform (initial and final stationary general eqms)
% of who becomes a worker vs entpreneur (in terms of assets and entrepreneurial ability; a and z)
% Note: You will need to make this figure full screen to be able to read it
% If we plot all the asset grid the ylabels is a mess because of the lines (between each square), so just plot every 10th asset grid point
% (better would plot all and then modify the y-axis so that doesn't ytick them all, but I can't be bothered)
if CreateFigures==1
    just10thasset=1:10:n_a;
    figure(3)
    temp1=gather(workerORentrepreneur_init(just10thasset,:,1,1)); % heatmap only works with cpu
    heatmap(z_grid,a_grid(just10thasset),temp1)
    grid off
    title('Initial eqm, with tax: Who becomes entrepreneur')
end
% Can see that for initial eqm with tax, only very high entrepreneurial ability become entrepreneurs, and only those with some assets

%% Save results to mat file
save BS2013_preTpath.mat

%% Replicate Figure 2 of BS2013

%ii_bench    = 1;
lambda_vec = [inf,2.0,1.75,1.5,1.25,1.0]';
NN = length(lambda_vec);
do_GE = 1;

%Pre-allocate arrays or structures where you want to store the output
share_entre_vec = zeros(NN,1);
extfin_vec      = zeros(NN,1);
extfin_Y_vec    = zeros(NN,1);
Y_vec           = zeros(NN,1);
r_vec           = zeros(NN,1);
w_vec           = zeros(NN,1);
K_vec           = zeros(NN,1);
K_Y_vec         = zeros(NN,1);

for ii=1:length(lambda_vec)

    %Assign lambda:
    Params.lambda = lambda_vec(ii);

    disp('***************************************************************')
    fprintf('Doing experiment %d of %d \n',ii,length(lambda_vec));
    fprintf('lambda = %.3f \n',lambda_vec(ii));
    disp('***************************************************************')

    if ii==1
        Params.r = 0.0472;
        Params.w = 0.171;
    elseif ii>1
        Params.r = r_vec(ii-1) ;
        Params.w = w_vec(ii-1);
    end

    tic
    [Outputs] = BueraShin_Fn(do_GE,Params,n_d,n_a,n_z,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
    toc

    %Aggregate quantities and prices
    share_entre_vec(ii) = Outputs.share_entre;
    extfin_vec(ii)      = Outputs.extfin;
    Y_vec(ii)           = Outputs.Y;
    r_vec(ii)           = Outputs.r;
    w_vec(ii)           = Outputs.w;
    K_vec(ii)           = Outputs.K;
    K_Y_vec(ii)         = Outputs.K_Y;
    extfin_Y_vec(ii)    = Outputs.extfin/Outputs.Y;

end

%Add "_norm" to denote change wrt benchmark
ii_bench    = 1;
Y_vec_norm  = zeros(NN,1);

for ii=1:NN
    Y_vec_norm(ii) = (Y_vec(ii)/Y_vec(ii_bench));
end

%% Plots for Figure 2 of the paper

ldw = 2;
fts = 14;

figure(1)
plot(extfin_Y_vec,Y_vec_norm,'linewidth',ldw)
xlabel('External Finance to GDP','FontSize',fts)
ylabel('GDP relative to benchmark','FontSize',fts)
title('GDP and TFP','FontSize',fts)
print('fig2a_BS2013','-dpng')

figure(2)
plot(extfin_Y_vec,r_vec,'linewidth',ldw)
xlabel('External Finance to GDP','FontSize',fts)
ylabel('Interest rate','FontSize',fts)
title('Interest Rate','FontSize',fts)
print('fig2b_BS2013','-dpng')

figure(3)
subplot(1,2,1)
    plot(extfin_Y_vec,Y_vec_norm,'linewidth',ldw)
    xlabel('External Finance to GDP','FontSize',fts)
    ylabel('GDP relative to benchmark','FontSize',fts)
    title('GDP and TFP','FontSize',fts)
subplot(1,2,2)
    plot(extfin_Y_vec,r_vec,'linewidth',ldw)
    xlabel('External Finance to GDP','FontSize',fts)
    ylabel('Interest rate','FontSize',fts)
    title('Interest Rate','FontSize',fts)
print('fig2_BS2013','-dpng')


save data_all 