%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pp] = mk_pointer_1(num)

ns = num.stress;
% nd = num.dof;

%%%% variables : [e]      \in (ns)                     %%%%
%%%% variables : [s]      \in (ns)                     %%%%
%%%% variables : [u]      \in (nd)                     %%%%

pp.e = 0;
pp.s = pp.e + ns;
pp.u = pp.s + ns;

