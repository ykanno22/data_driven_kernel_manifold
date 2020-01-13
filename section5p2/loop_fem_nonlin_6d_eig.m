%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
% --->
global param num vec_a cons_c data_x lin_eq
load('nonlin_eig_6d.mat');
% <---
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh size, etc.
% --->
num_Iter = 10;
%
Length = 60.0; % Length of the model in [mm]
Width  = 20.0; % Width               in [mm]
NXE = 18;       % Number of rows in the x direction
NYE = 6;       % Number of rows in the y direction
%
nd = 2 * NXE * (NYE+1);
idx_bottom_right_ydirec = nd - (2 * NYE);
%
ngp  = 2;
% <---
% Mesh size, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
his_load_fact = zeros(1, num_Iter+1);
his_u_mean = zeros(nd, num_Iter+1);
his_sol_u  = zeros(nd, num_Iter+1);

fprintf(' ============================================= \n');
fprintf(' size of elastic body = %3.2f mm x %3.2f mm \n',...
    Length, Width);
fprintf(' for conventional FEM:   Young = %3.5f MPa ;  Poisson = %3.5f \n',...
    mean_Young * param_scale, mean_Poisson);
fprintf(' compare displacement at right bottom node, y-direction \n');

incr_load_factor = 1 / num_Iter;
sol_u   = zeros(nd, 1);
sol_eps = zeros(3 * NXE * NYE * ngp * ngp, 1);
sol_sig = zeros(3 * NXE * NYE * ngp * ngp, 1);

fprintf('     ----------------------------------------- \n');
for Iter = 0:num_Iter
    load_fact = Iter * incr_load_factor;
    
    [vec_u_mean, sol_u, sol_eps, sol_sig, output] =...
        fun_fem_nonlin_6d_eig(NXE, NYE, Length, Width,...
        sol_u, sol_eps, sol_sig, load_fact);
    
    his_u_mean(:,Iter+1) = vec_u_mean;
    his_sol_u(:,Iter+1)  = sol_u;
    his_load_fact(1,Iter+1) = load_fact;
    prev_u = sol_u;
    
    fprintf('         #iter = %g ;  #func_call = %g \n',...
        output.iterations, output.funcCount );
    fprintf('     ----------------------------------------- \n');
end

figure;
plot(his_u_mean(idx_bottom_right_ydirec,:), his_load_fact, 'bo:',...
    'LineWidth', 1.0, 'MarkerSize',13);
hold on;
plot(his_sol_u(idx_bottom_right_ydirec,:), his_load_fact, 'rs-',...
    'LineWidth', 1.0, 'MarkerSize',10);
xlabel('Displacement (mm)', 'Interpreter', 'latex');
ylabel('Load factor, $\lambda$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);

