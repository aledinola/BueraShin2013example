function capitaldemand=BueraShin2013_capitaldemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

[profit,kstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

if w>profit
    capitaldemand=0; % worker
else
    capitaldemand=kstar; % entrepreneur
end

end %end function