function [exit_E_to_W,entry_W_to_E,T] = fun_entry_exit(StationaryDist,pi_z,pol_e,pol_aprime,n_a,n_z)
% This small function computes the entry and exit rates b/w occupations. 
% INPUTS:
% StationaryDist: 
%   pol_e:      Policy function for entre choice (0=worker,1=entre)
%   pol_aprime: Policy function for assets, indexes
%   pi_z:       Transition matrix prob(z,z')

% Check inputs
if ~isequal(size(StationaryDist),[n_a,n_z])
    error('Size of StationaryDist NOT correct!')
end
if ~isequal(size(pol_e),[n_a,n_z])
    error('Size of pol_e NOT correct!')
end
if ~isequal(size(pol_aprime),[n_a,n_z])
    error('Size of pol_aprime NOT correct!')
end
if ~isequal(size(pi_z),[n_z,n_z])
    error('Size of pi_z NOT correct!')
end 

% Compute T(o,o')=mass of individuals who go from occupation o in the current 
% period to occupation o' in the next period

T = zeros(2,2); %1=worker, 2=entre
 
for z_c=1:n_z
    for a_c=1:n_a
        aprime_c = pol_aprime(a_c,z_c);
        occ = pol_e(a_c,z_c);
        o = 1*(occ==0)+2*(occ==1); %1=worker, 2=entre
        for zprime_c = 1:n_z
            occ_prime = pol_e(aprime_c,zprime_c);
            o_prime   = 1*(occ_prime==0)+2*(occ_prime==1); %1=worker, 2=entre
            T(o,o_prime) = T(o,o_prime)+pi_z(z_c,zprime_c)*StationaryDist(a_c,z_c);
        end
    end
end

exit_E_to_W = T(2,1)/sum(T(2,:));
entry_W_to_E = T(1,2)/sum(T(1,:));

end %end function