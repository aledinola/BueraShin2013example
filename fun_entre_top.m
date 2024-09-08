function [top10_empl] = fun_entre_top(pol_labor,pol_e,StationaryDist,n_a,n_z)
%This function computes one of the moments reported in Table 1 of BS2013:
% (1) Top 10% of employment

% Check inputs
if ~isequal(size(StationaryDist),[n_a,n_z])
    error('Size of StationaryDist NOT correct!')
end
if ~isequal(size(pol_labor),[n_a,n_z])
    error('Size of pol_labor NOT correct!')
end
if ~isequal(size(pol_e),[n_a,n_z])
    error('Size of pol_e NOT correct!')
end

% Consider entrepreneurs only. If (a,z) is a worker, than labor=0
lab    = pol_labor(:);      %(n_a*n_z,1)
entre  = pol_e(:);          %(n_a*n_z,1)
mass   = StationaryDist(:); %(n_a*n_z,1)
lab(entre==0)  = [];
mass(entre==0) = [];
mass = mass/sum(mass);

[lab_sorted,ind_sorted] = sort(lab);
mass_sorted = mass(ind_sorted);
cdf_labor   = cumsum(mass_sorted);
[~,p_90]    = min(abs(cdf_labor-0.90));
numerator   = sum(lab_sorted(p_90+1:end).*mass_sorted(p_90+1:end));
denominator = sum(lab_sorted.*mass_sorted);
top10_empl  = numerator/denominator;

end %end function