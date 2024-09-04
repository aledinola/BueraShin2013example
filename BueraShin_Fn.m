function [Outputs,Policy,StationaryDist,AggVars,ValuesOnGrid] = BueraShin_Fn(Params,n_d,n_a,n_z,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% Note: The initial guesses for the equilibrium prices r and w are stored
% in Params.r and Params.w

%% Compute GE
[p_eqm,~,~]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

%% Compute eqm objects
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, [], simoptions);

%% Compute some model moments

%Aggregate quantities and prices
Outputs.share_entre = AggVars.entrepreneur.Mean;
Outputs.extfin      = AggVars.extfin.Mean;
Outputs.Y           = AggVars.Y.Mean;
Outputs.r           = Params.r;
Outputs.w           = Params.w;
Outputs.K           = AggVars.K.Mean;
Outputs.K_Y         = Outputs.K/Outputs.Y;
Outputs.extfin_Y    = Outputs.extfin/Outputs.Y;

end %end function BueraShin_Fn