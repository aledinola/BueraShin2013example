function entrepreneur=BueraShin2013_entrepreneur(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

profit = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

if w>profit
    entrepreneur=0;
else
    entrepreneur=1;
end

end %end function