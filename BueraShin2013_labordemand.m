function labordemand=BueraShin2013_labordemand(aprime,a,z,w,r,lambda,delta,alpha,upsilon)

[profit,~,lstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);


if w>profit
    labordemand=0; % worker
else
    labordemand=lstar; % entrepreneur
end

end %end function