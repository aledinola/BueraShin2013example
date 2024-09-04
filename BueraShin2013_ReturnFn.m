function F=BueraShin2013_ReturnFn(aprime,a,z,crra,w,r,lambda,delta,alpha,upsilon)

F=-Inf;

profit = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

% Budget constraint
c=max(w,profit)+(1+r)*a-aprime;

if c>0
    F=(c^(1-crra))/(1-crra);
end

end %end function