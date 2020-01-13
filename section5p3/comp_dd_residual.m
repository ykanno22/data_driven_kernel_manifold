%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = comp_dd_residual(x)

global param num vec_a cons_c data_x data_cur lin_eq

ns = num.stress;
np = num.plane;

F = [lin_eq.A * x - lin_eq.b];

vec_eps = x(1:ns);
vec_sig = x(ns+1:(2*ns));
pp = 0;
on_plane = zeros(ns,1);

for i=1:ns/3
    eps_i = vec_eps(pp+1:pp+3);
    sig_i = vec_sig(pp+1:pp+3);
    x_i = [eps_i; sig_i];
    for k=1:np
        vec_a_at_x  = zeros(num.dim,1);
        cons_c_at_x = 0;
        for l=1:num.data
            x_l = data_x(1:3,l);
            kappa_l = exp(-param.RBF * ((eps_i-x_l)'*(eps_i-x_l)) );
            vec_a_at_x  = vec_a_at_x  + (kappa_l * vec_a{k,l});
            cons_c_at_x = cons_c_at_x + (kappa_l * cons_c{k,l});
        end
        pp = pp + 1;
        on_plane(pp) = (vec_a_at_x' * x_i) - cons_c_at_x;
    end
end

F = [F; on_plane];

