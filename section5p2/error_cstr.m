function [c, ceq] = error_cstr(x)

global param num vec_a cons_c data_x

d_x = num.dim_x;
on_plane = zeros(param.error_cur,1);
for k=1:param.error_cur
    vec_a_at_x  = zeros(num.dim,1);
    cons_c_at_x = 0;
    for l=1:num.data
        x_l = data_x(:,l);
        kappa_l = exp(-param.RBF * ((x(1:d_x)-x_l(1:d_x))'*(x(1:d_x)-x_l(1:d_x))) );
        vec_a_at_x  = vec_a_at_x  + (kappa_l * vec_a{k,l});
        cons_c_at_x = cons_c_at_x + (kappa_l * cons_c{k,l});
    end
    on_plane(k) = (vec_a_at_x' * x) - cons_c_at_x;
end

ceq = on_plane;
c = [];