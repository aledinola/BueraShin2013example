clear
clc
close all
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))
% Buera & Shin (2013) - Financial Frictions and the Persistence of History: 
% A Quantitative Exploration
% The model of Buera & Shin (2013) requires some reworking before
% computing, see the accompanying pdf.
% I use z to denote the markov shock (as that is standard nomenclature of
% VFI Toolkit), which was called e by BS2013.

%% Set folders
ResFolder = 'results'; % Folder to save results

%% Flags and numerical options
CreateFigures  = 1; % Flag 0/1 plot figures of initial steady-state
do_GE          = 1; % 0 = partial equilibrium, 1 = general equilibrium
do_replication = 0; % Flag 0/1 to replicate Figure 2 of BS 2013. This 
                    % requires repeatedly solving the s.s. for different
                    % lambdas
heteroagentoptions.fminalgo=1;
heteroagentoptions.maxiter=50;
heteroagentoptions.verbose=1;
heteroagentoptions.toleranceGEprices=10^(-3);
heteroagentoptions.toleranceGEcondns=10^(-3);

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
Params.lambda  =inf; % Calibration of steady-state done for US economy

% Initial values for general eqm parameters: good for lambda=inf
Params.r = 0.0476;
Params.w = 0.172;

%% Grid for assets
d_grid = []; % No grid for static choice d
n_d    = 0;  % No grid for static choice d
% grid on assets
n_a    = 1001; % Num of grid points
a_min  = 1e-6; % Lower bound
a_max  = 4000; % Upper bound
% a_scale>1 puts more points near zero
a_scale = 2;   % "Curvature" of asset grid
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
n_z      = length(z_grid); % Num of grid points for exo state z
% See my notes on Buera and Shin paper
pi_z = Params.psi*eye(n_z)+(1-Params.psi)*ones(n_z,1)*pi_z_vec';
pi_z = pi_z./sum(pi_z,2);

%% Return fn
DiscountFactorParamNames={'beta'};
% Required inputs:
% (aprime,a,z) in this order, than any parameter
ReturnFn=@(aprime,a,z,crra,w,r,lambda,delta,alpha,upsilon) ...
    BueraShin2013_ReturnFn(aprime,a,z,crra,w,r,lambda,delta,alpha,upsilon);

%% Create some FnsToEvaluate
FnsToEvaluate.A=@(aprime,a,z) a; % assets
% Capital demand by entre (zero if worker)
FnsToEvaluate.K=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_capitaldemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% Labor demand by entre (zero if worker)
FnsToEvaluate.L=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_labordemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% 1 if entrepreneur, 0 if worker
FnsToEvaluate.entrepreneur=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_entrepreneur(aprime,a,z,w,r,lambda,delta,alpha,upsilon);
% Entrepreneurial output (zero if worker)
FnsToEvaluate.Y=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_output(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% Enternal finance (zero if worker)
FnsToEvaluate.extfin=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_extfin(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 
% Earnings 
FnsToEvaluate.earnings=@(aprime,a,z,w,r,lambda,delta,alpha,upsilon) BueraShin2013_earnings(aprime,a,z,w,r,lambda,delta,alpha,upsilon); 

%% Model is set, we can start by just that the basics are running okay
vfoptions           = struct();
vfoptions.verbose   = 0;
vfoptions.lowmemory = 0;
simoptions          = struct();

%% Set up general equilibrium
% heteroagentoptions.fminalgo=5;
% % Need to explain to heteroagentoptions how to use the GeneralEqmEqns to update the general eqm prices.
% heteroagentoptions.fminalgo5.howtoupdate={...  % a row is: GEcondn, price, add, factor
%     'capitalmarket','r',1,0.03;...  % capitalmarket GE condition will be positive if r is too big, so subtract
%     'labormarket','w',1,0.1;... % labormarket GE condition will be positive if w is too small, so add
% };

% There are two prices to be determined in GE and two GE conditions
GEPriceParamNames={'r','w'};
% GE(1): capital demand minus capital supply
GeneralEqmEqns.capitalmarket=@(K,A) K-A; 
% GE(2): labor demand minus labor supply, suppy is just fraction of workers (who each exogneously supply endowment 1 of labor)
GeneralEqmEqns.labormarket=@(L,entrepreneur) L-(1-entrepreneur); 

%% Compute the model once, either in partial or general equilibrium

if do_GE==1
    disp('Solving initial stationary general eqm')
else
    disp('Solving initial partial eqm')
end

tic
[Outputs,GE_cond,Policy_init,StationaryDist_init,AggVars_init,ValuesOnGrid] = BueraShin_Fn(do_GE,Params,n_d,n_a,n_z,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
time_BueraShin_Fn=toc;

disp(' ')
fprintf('time_BueraShin_Fn = %f\n',time_BueraShin_Fn)

%% Analyse results from model solution and make plots
Params.r=Outputs.r;
Params.w=Outputs.w;

workerORentrepreneur_init=ValuesOnGrid.entrepreneur;
clear ValuesOnGrid

if CreateFigures==1
    % take another a look at cdf over asset grid to make sure not hitting top of grid
    figure
    plot(a_grid,cumsum(sum(sum(sum(StationaryDist_init,4),3),2)))
    title('cdf of asset to make sure grid on assets seems okay (init eqm)')

    just10thasset=1:10:n_a;
    figure
    temp1=gather(workerORentrepreneur_init(just10thasset,:,1,1)); % heatmap only works with cpu
    heatmap(z_grid,a_grid(just10thasset),temp1)
    grid off
    title('Initial eqm, with tax: Who becomes entrepreneur')
end

%% Replicate Table 1 of BS2013

make_table(ResFolder,GE_cond,Outputs);

%% Replicate Figure 2 of BS2013

if do_replication==1
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

figure
plot(extfin_Y_vec,Y_vec_norm,'linewidth',ldw)
xlabel('External Finance to GDP','FontSize',fts)
ylabel('GDP relative to benchmark','FontSize',fts)
title('GDP and TFP','FontSize',fts)
print('fig2a_BS2013','-dpng')

figure
plot(extfin_Y_vec,r_vec,'linewidth',ldw)
xlabel('External Finance to GDP','FontSize',fts)
ylabel('Interest rate','FontSize',fts)
title('Interest Rate','FontSize',fts)
print(fullfile(ResFolder,'fig2b_BS2013'),'-dpng')

figure
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
print(fullfile(ResFolder,'fig2_BS2013'),'-dpng')

save (fullfile(ResFolder,"data_all.mat")) 

end %end if do_replication