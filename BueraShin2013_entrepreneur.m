function entrepreneur=BueraShin2013_entrepreneur(aprime,a,z,tau,psi,w,r,lambda,delta,alpha,upsilon)

% Get k1, kstar, lstar
k1a=(1/(r+delta))*alpha*(1-upsilon)*(1-tau)*z;
k1b=(1/w)*(1-alpha)*(1-upsilon)*(1-tau)*z;
k1=(k1a^(1-(1-alpha)*(1-upsilon)) * k1b^((1-alpha)*(1-upsilon)))^(1/upsilon);
kstar=min(k1,lambda*a);
lstar=( (1/w)*(1-alpha)*(1-upsilon)*(1-tau)*z *kstar^(alpha*(1-upsilon)) )^(1/(1-(1-alpha)*(1-upsilon)));
% Evaluate profit if do choose to be entrepreneur
pi=(1-tau)*z*((kstar^alpha)*(lstar^(1-alpha)) )^(1-upsilon) -w*lstar -(delta+r)*kstar;

if w>pi
    entrepreneur=0;
else
    entrepreneur=1;
end


end