function output=BueraShin2013_output(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

[profit,kstar,lstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

if w>profit
    output=0; % worker
else
    output=z*((kstar^alpha)*(lstar^(1-alpha)) )^(1-upsilon); % entrepreneur
end

end %end function