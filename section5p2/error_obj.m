function [f_val, grad_f] = error_obj(x)

global data_cur

f_val = (x - data_cur)' * (x - data_cur);

if nargout > 1
    grad_f = 2 * (x - data_cur);
end