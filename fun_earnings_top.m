function [top5_earnings] = fun_earnings_top(pol_earnings,StationaryDist,n_a,n_z)
%This function computes one of the moments reported in Table 1 of BS2013:
% (1) top 5% earnings

% Check inputs
if ~isequal(size(StationaryDist),[n_a,n_z])
    error('Size of StationaryDist NOT correct!')
end
if ~isequal(size(pol_earnings),[n_a,n_z])
    error('Size of pol_earnings NOT correct!')
end

% Vectorize
value = pol_earnings(:);   %(n_a*n_z,1)
mass  = StationaryDist(:); %(n_a*n_z,1)

[value_sorted,ind_sorted] = sort(value);
mass_sorted = mass(ind_sorted);
cdf         = cumsum(mass_sorted);
[~,p_95]    = min(abs(cdf-0.95));
numerator   = sum(value_sorted(p_95+1:end).*mass_sorted(p_95+1:end));
denominator = sum(value_sorted.*mass_sorted);
top5_earnings = numerator/denominator;

end %end function