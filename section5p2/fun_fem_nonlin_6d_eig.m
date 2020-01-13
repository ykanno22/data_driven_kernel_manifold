%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%% This code is partially based on:                         %%%%%%%%%%
%%%%%   A. Khennane:                                           %%%%%%%%%%
%%%%%   Introduction to Finite Element Analysis Using MATLAB and Abaqus.%
%%%%%   CRC Press, Boca Raton (2013).                          %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_u_mean, sol_u, sol_eps, sol_sig, output] =...
    fun_fem_nonlin_6d_eig(NXE, NYE, Length, Width,...
    prev_u, prev_eps, prev_sigma, load_fact)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
% --->
global param num vec_a cons_c data_x lin_eq
load('nonlin_eig_6d.mat');
% <---
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To change the size of the mesh, alter the next statements
% %%%%% Introduction to Finite Element Analysis Using MATLAB and Abaqus
% --->
thick = 5.0;      % Tickness in mm
% E  = 20.0;      % Elastic modulus in  100*MPa
% nu = 0.3;      % Poisson's ratio 
E  = mean_Young;
nu = mean_Poisson;
%
% Force = 200.0;    % external force in kN
Force = load_fact * (0.25 / NXE); %%% distributing: Force / (Length * thhick) in [kN/mm2]
%
dhx = Length/NXE; % Element size in the x direction
dhy = Width/NYE;  % Element size in the x direction
X_origin = 0.0;       % X origin of the global coordinate system
Y_origin = 0.0;   % Y origin of the global coordinate system
% <--
% To change the size of the mesh, alter the next statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the mesh 
% --->
nne = 4;
nodof = 2;
eldof = nne*nodof;
%
[nnd, nel, geom, connec] =Q4_mesh_fun(NXE, NYE, X_origin, Y_origin, dhx, dhy); 
% <--
% Generate the mesh 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the elastic matrix for plane stress 
% --->
dee = formdsig(E, nu);
% <--
% Form the elastic matrix for plane stress 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions
% --->
%%%% Restrain in all directions the nodes situated @ (x = 0)
nf = ones(nnd, nodof);
%
fixed_dig = (geom(:,1)==0);
nf(fixed_dig,:) = 0;
clear fixed_dig
%
nd = sum(sum(nf));
idx_bottom_right_ydirec = nd - (2 * NYE);
%
ii = 0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0
            ii = ii + 1;
            nf(i,j) = ii;
        end
    end
end

%%%% Apply uniform distributed load at the nodes at y=Width
Nodal_loads= zeros(nnd, 2);
Ind_top_nodes = find(geom(:,2) == Width);
Nodal_loads(Ind_top_nodes,2)      = -Force;
Nodal_loads(Ind_top_nodes(end),2) = -Force/2;
%
Nodal_loads = sparse(Nodal_loads);
% <--
% Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the global force vector
% --->
vec_f = zeros(nd,1);
%
Idx_tmp = find(nf(:,1));
vec_f(nf(Idx_tmp,1)) = Nodal_loads(Idx_tmp,1);
clear Idx_tmp
%
Idx_tmp = find(nf(:,2));
vec_f(nf(Idx_tmp,2)) = Nodal_loads(Idx_tmp,2);
clear Idx_tmp
%
vec_f = sparse(vec_f) * param_scale;
% <---
% Assemble the global force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the matrix containing the abscissas and the weights of Gauss points
% --->
ngp  = 2;
samp = gauss(ngp);
% <---
% Form the matrix containing the abscissas and the weights of Gauss points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration and assembly of the global stiffness matrix
% --->
mat_B = [];
integ_weight = [];
for i=1:nel
    %%%% coordinates of the nodes of element i, & its steering vector -->
    [coord, g] = elem_q4_fun(i, nne, nodof, geom, connec, nf);
    ke = sparse(eldof, eldof);
    for ig=1:ngp
        wi = samp(ig,2);
        for jg=1:ngp
            wj = samp(jg, 2);
            %%%% Derivative of shape functions in local coordinates -->
            [der, fun] = fmlin(samp, ig, jg);
            %%%% Compute Jacobian matrix -->
            jac = der * coord;
            %%%% Compute determinant of Jacobian matrix -->
            d = det(jac);
            %%%% Derivative of shape functions in global coordinates -->
            deriv = jac \ der;
            %%%% Form matrix [B] -->
            bee = formbee(deriv, nne, eldof);
            %%%% Integrate stiffness matrix -->
            integ_weight = [integ_weight, d * thick * wi * wj];
            mat_Be = zeros(3, nd);
            for jj=find(g)
                mat_Be(:,g(jj)) = bee(:,jj);
            end
            mat_B = [mat_B; mat_Be];
        end
    end
end
mat_K = mat_B' * kron(diag(integ_weight), dee) * mat_B;
% <---
% Numerical integration and assembly of the global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution process
% --->
vec_u_mean = mat_K \ vec_f;
% <---
% Solution process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% variables : [e]      \in (ns)                     %%%%
%%%% variables : [s]      \in (ns)                     %%%%
%%%% variables : [u]      \in (nd)                     %%%%
num.dof = nd;
nm = nel;     num.element = nm;
ns = 3 * (ngp*ngp) * nm;  num.stress = ns;  num.strain = ns;
num.all_var = ns + ns + nd;

lin_eq.A = [];
lin_eq.b = [];

%%% equalities for (s) : force-balance eq.
pp = mk_pointer_1(num);

Aeq_f = zeros(nd,num.all_var);
Aeq_f(:,(pp.s+1):(pp.s+ns)) =  (kron(diag(integ_weight), eye(3)) * mat_B)';
Aeq_f = sparse(Aeq_f);
beq_f = vec_f;
%
lin_eq.A = [lin_eq.A; Aeq_f];
lin_eq.b = [lin_eq.b; beq_f];
clear Aeq_f beq_f

%%% equalities for (e, u) : compatibility eq.
pp = mk_pointer_1(num);

Aeq_c = zeros(ns,num.all_var);
Aeq_c(:,(pp.u+1):(pp.u+nd)) = -mat_B;
Aeq_c(:,(pp.e+1):(pp.e+ns)) =  speye(ns);
Aeq_c = sparse(Aeq_c);
beq_c = sparse(ns,1);
%
lin_eq.A = [lin_eq.A; Aeq_c];
lin_eq.b = [lin_eq.b; beq_c];
clear Aeq_c beq_c



x0 = [prev_eps; prev_sigma; prev_u];


% options = optimoptions('fsolve', 'MaxIter',2000,...
%     'Algorithm','trust-region-dogleg', 'FinDiffType','central',...
%     'Display','off');
options = optimoptions('fsolve', 'MaxIter',2000, 'MaxFunEvals',10^8,...
    'Display','off');

[x,~,exitflag, output] = fsolve(@comp_dd_residual, x0, options);

[sol_eps, sol_sig, sol_u] = mk_out_1(x, num);



fprintf('     load factor = %1.5f \n',...
    load_fact );
fprintf('     conventional FEM : %3.5f mm\n',...
    vec_u_mean(idx_bottom_right_ydirec) );
fprintf('     data-driven solver : %3.5f mm  w/ exitflag = %g \n',...
    sol_u(idx_bottom_right_ydirec), exitflag);


