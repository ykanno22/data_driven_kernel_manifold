%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = comp_model_residual(x)

global param num data_x lin_eq mean_Young mean_Poisson

ns = num.stress;

F = [lin_eq.A * x - lin_eq.b];

vec_eps = x(1:ns);
vec_sig = x(ns+1:(2*ns));
pp = 0;
resid_material = zeros(ns,1);

for i=1:ns/3
    eps_i = vec_eps(pp+1:pp+3);
    sig_i = vec_sig(pp+1:pp+3);
    %
    mat_eps_i = [...
        eps_i(1),   eps_i(3)/2, 0;...
        eps_i(3)/2, eps_i(2),   0;...
        0, 0, (-mean_Poisson/(1-mean_Poisson)) * (eps_i(1) + eps_i(2))];
    dev_eps_i = mat_eps_i - ( (trace(mat_eps_i) / 3) * eye(3) );
    norm_d_eps_i = norm(dev_eps_i, 'fro');
    %
    Young_i = mean_Young...
        * ((3/2) - (1 / (1+exp(-10*norm_d_eps_i))));
    dee = (  Young_i / (1 - (mean_Poisson^2))  ) *...
        [1, mean_Poisson, 0;...
        mean_Poisson, 1, 0;...
        0, 0, (1/2) * (1 - mean_Poisson)];
    resid_material(pp+1:pp+3) = sig_i - (dee * eps_i);
    pp = pp + 3;
end

F = [F; resid_material];

