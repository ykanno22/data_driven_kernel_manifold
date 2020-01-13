function [vec_eps, vec_sig, vec_u] = mk_out_1(x, num)

ns = num.stress;
nd = num.dof;

%%%% variables : [e]      \in (ns)                     %%%%
%%%% variables : [s]      \in (ns)                     %%%%
%%%% variables : [u]      \in (nd)                     %%%%

vec_eps = x(1:ns);
x = x(ns+1:end);

vec_sig = x(1:ns);
x = x(ns+1:end);

vec_u = x(1:nd);
x = x(nd+1:end);

if ~isempty(x)
    fprintf('\n !!!!! Error in mk_out_1.m !!!!! \n');
end

