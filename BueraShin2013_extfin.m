function extfin=BueraShin2013_extfin(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

[profit,kstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);


if w>profit
    extfin=0; % worker
else
    extfin=max(0,kstar-a); % entrepreneur
end

end %end function