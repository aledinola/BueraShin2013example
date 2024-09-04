clear,clc,close all

% set params
gamma   = 1.5; % CES utility param
delta   = 0.06;
alpha   = 0.33;
upsilon = 0.21; % span of control = 1-upsilon
lambda  = 1.35;
% set prices
r = 0.04;
w = 1.1;
% set assets and ability
a = 18;
z = 1.3;

% set grids for k and l
n_k = 3000;
n_l = 1000;
k_grid = linspace(0,80,n_k)';
l_grid = linspace(0,50,n_l)';

% Compute profit(k,l) on the grid --> profit_mat has size [n_k,nl]
profit_mat = z*((k_grid.^alpha).*(l_grid'.^(1-alpha))).^(1-upsilon)-w*l_grid'-(delta+r)*k_grid;
not_allowed = repmat(k_grid>lambda*a,1,n_l);
profit_mat(not_allowed) = -inf;

% Maximization wrt k and l on the grids
[max_val,max_ind] = max(profit_mat(:));
prof1 = max_val;
[k1,l1] = ind2sub(size(profit_mat),max_ind);
k1 = k_grid(k1);
l1 = l_grid(l1);

disp('MAXIM ON GRID')
fprintf('Optimal profit  = %f \n',prof1)
fprintf('Optimal capital = %f \n',k1)
fprintf('Optimal labor   = %f \n',l1)

% Use analytical results
[profit,kstar,lstar] = solve_entre(a,z,w,r,lambda,delta,alpha,upsilon);

disp('MAXIM ANALYTICAL ')
fprintf('Optimal profit  = %f \n',profit)
fprintf('Optimal capital = %f \n',kstar)
fprintf('Optimal labor   = %f \n',lstar)