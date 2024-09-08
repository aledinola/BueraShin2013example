function earnings=BueraShin2013_earnings(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

profit = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

% Earnings
earnings = max(w,profit) + r*a;


end %end function