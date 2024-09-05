function [profit,kstar,lstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon)
% This function solve the entrepreneurs' static maximization problem. See
% Robert's App_BueraShin2013.pdf for more details.

% Get k1, kstar, lstar
%aux    = 1-(1-alpha)*(1-upsilon);
k1a    = (1/(max(r+delta,1e-8)))*alpha*(1-upsilon)*z;
k1b    = (1/w)*(1-alpha)*(1-upsilon)*z;
inside = k1a^(1-(1-alpha)*(1-upsilon)) * k1b^((1-alpha)*(1-upsilon));
k1     = (inside)^(1/upsilon);
kstar  = min(k1,lambda*a);
inside_lab = (1/w)*(1-alpha)*(1-upsilon)*z *kstar^(alpha*(1-upsilon));
lstar  = ( inside_lab )^(1/(1-(1-alpha)*(1-upsilon)));
% Evaluate profit if do choose to be entrepreneur
profit=z*((kstar^alpha)*(lstar^(1-alpha)) )^(1-upsilon) -w*lstar -(delta+r)*kstar;

end %end function