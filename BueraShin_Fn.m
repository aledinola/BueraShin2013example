function [Outputs,GE_cond,Policy,StationaryDist,AggVars,ValuesOnGrid] = BueraShin_Fn(do_GE,Params,n_d,n_a,n_z,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% Note: The initial guesses for the equilibrium prices r and w are stored
% in Params.r and Params.w

if do_GE==1
    % When doing GE, better not to display intermediate output from VFI
    vfoptions.verbose=0;
end

if do_GE==1
    %% Compute GE
    [p_eqm,~,~]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

    %% Compute eqm objects
    Params.r=p_eqm.r;
    Params.w=p_eqm.w;
    [~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
else
    [~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
    ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
end

%% Compute GE conditions

% Excess demand capital market: capital demand minus capital supply
GE_cond(1) = AggVars.K.Mean-AggVars.A.Mean;
% Excess demand labor market: labor demand minus share of workers
GE_cond(2) = AggVars.L.Mean - (1-AggVars.entrepreneur.Mean);

%% Compute some model moments using the toolkit options
% NOTE: I prefer to do my own moment calculations since the toolkit
% commands are a bit slow. I verified that I get almost identical results

% % - Conditional restrictions
% % Condition=1 if ENTRE, 0 if WORKER
% simoptions.conditionalrestrictions.ENT = FnsToEvaluate.entrepreneur;
% 
% simoptions.whichstats=zeros(1,7);
% simoptions.whichstats(1)=0; %mean
% simoptions.whichstats(2)=0; %median
% simoptions.whichstats(4)=1; %lorenz curve and gini
% simoptions.whichstats(6)=0; %quantiles (1=loops,2=fast but lots of RAM)
% %simoptions.nquantiles = 100;
% 
% AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate,...
%     Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
% 
% % --- Earnings distribution: top shares
% %top1_earnings  = 1-AllStats.A.LorenzCurve(99);
% top5_earnings2  = 1-AllStats.earnings.LorenzCurve(95);
% %top10_earnings = 1-AllStats.A.LorenzCurve(90);
% 
% % --- Employment distribution (conditional on entre status): top shares
% top10_empl2 = 1-AllStats.ENT.L.LorenzCurve(90);

%% Compute some model moments using my own calculations

pol_e      = gather(ValuesOnGrid.entrepreneur); % dim: (a,z)
pol_aprime = gather(squeeze(Policy(1,:,:))); % dim: (a,z)
pol_labor  = gather(ValuesOnGrid.L); % dim: (a,z)
pol_earnings = gather(ValuesOnGrid.earnings); % dim: (a,z)

% Compute exit rate
[exit_E_to_W,entry_W_to_E] = fun_entry_exit(StationaryDist,pi_z,pol_e,pol_aprime,n_a,n_z);

% Compute top 10% employment (i.e. share of employment in the largest 10%
% of establishments)
[top10_empl] = fun_entre_top(pol_labor,pol_e,StationaryDist,n_a,n_z);

% Compute top 5% earnings
[top5_earnings] = fun_earnings_top(pol_earnings,StationaryDist,n_a,n_z);

%Aggregate quantities and prices
Outputs.share_entre = AggVars.entrepreneur.Mean;
Outputs.extfin      = AggVars.extfin.Mean;
Outputs.Y           = AggVars.Y.Mean;
Outputs.r           = Params.r;
Outputs.w           = Params.w;
Outputs.K           = AggVars.K.Mean;
Outputs.K_Y         = Outputs.K/Outputs.Y;
Outputs.extfin_Y    = Outputs.extfin/Outputs.Y;
Outputs.exit_E_to_W   = exit_E_to_W;
Outputs.entry_W_to_E  = entry_W_to_E;
Outputs.top10_empl    = top10_empl;
Outputs.top5_earnings = top5_earnings;

% % Debug
% fprintf('top10_empl = %f \n',top10_empl)
% fprintf('top10_empl2 = %f \n',top10_empl2)
% fprintf('top5_earnings = %f \n',top5_earnings)
% fprintf('top5_earnings2 = %f \n',top5_earnings2)

end %end function BueraShin_Fn