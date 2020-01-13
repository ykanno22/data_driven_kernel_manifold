%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
% --->
global param num vec_a cons_c data_x data_cur
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.RBF   = 0.05 * 10^(0);
param.regularize = 10.0 * 10^(1);
%
load('3d_line_data.mat');
num.data = length(list_of_noisy_x);
num.dim  = 3;
num.variable = num.data * (num.dim + 1);
%
data_x(1,:) = list_of_noisy_x';
data_x(2,:) = list_of_noisy_y';
data_x(3,:) = list_of_noisy_f';
clear list_of_noisy_x list_of_noisy_y list_of_noisy_f
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare matrices & vectors
% - -->
matKappa = zeros(num.data,num.data);
vec_h = cell(1,num.data);
for i=1:num.data
    x_i = data_x(:,i);
    vec_h{i} = [];
    for j=1:num.data
        x_j = data_x(:,j);
        kappa_ij = exp(-param.RBF * ((x_i(1:2)-x_j(1:2))'*(x_i(1:2)-x_j(1:2))) );
        vec_h{i} = [vec_h{i};...
            kappa_ij * x_i; -kappa_ij];
        matKappa(i,j) = kappa_ij;
    end
end
matH  = zeros(num.data, num.variable);
for i=1:num.data
    matH(i,:) = vec_h{i}';
end

matR  = zeros(num.variable, num.variable);
for i=1:num.data
    x_i = data_x(:,i);
    for j=(i+1):num.data
        x_j = data_x(:,j);
        kappa_ij = exp(-1.0*param.RBF * ((x_i(1:2)-x_j(1:2))'*(x_i(1:2)-x_j(1:2))) );
        %
        pp.i = (i-1) * (num.dim+1);
        pp.j = (j-1) * (num.dim+1);
        d_ij = zeros(num.variable, num.dim+1);
        d_ij(pp.i + (1:(num.dim+1)), :) =  eye(num.dim+1);
        d_ij(pp.j + (1:(num.dim+1)), :) = -eye(num.dim+1);
        d_ij = sparse(d_ij);
        matR = matR + (kappa_ij * d_ij * (d_ij'));
        clear pp
    end
end

matQ = (matH' * matH) + (param.regularize * matR);
% <---
% Prepare matrices & vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning via eigenvalue analysis
% --->
[eig_vec, eig_val] = eig(matQ);
eig_val = diag(eig_val);
[eig_val, idx_sort] = sort(eig_val);
eig_vec = eig_vec(:,idx_sort);
%
fprintf(' ============================================= \n');
fprintf('   eigenvalue     square.error   regularization \n');
for i=1:7
    fprintf('     %3.5d    %3.5d    %3.5d \n',...
        eig_val(i),...
        eig_vec(:,i)' * (matH') * matH * eig_vec(:,i),...
        eig_vec(:,i)' * matR * eig_vec(:,i));
end
fprintf(' ============================================= \n');
%
vec_sol = cell(1,3);
for k=1:3
    vec_sol{k} = eig_vec(:,k);
end
%
vec_a  = cell(3,num.data);
cons_c = cell(3,num.data);
for k=1:3
    pp.cur = 0;
    for l = 1:num.data
        vec_a{k,l}  = vec_sol{k}(pp.cur + (1:num.dim));
        cons_c{k,l} = vec_sol{k}(pp.cur + num.dim+1);
        pp.cur = pp.cur + (num.dim + 1);
    end
end
% <---
% Learning via eigenvalue analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare figures: plane
% --->
x_plot = (min(data_x(1,:)) : max(data_x(1,:))/40 :max(data_x(1,:)));
x_plot = 1.05*[x_plot, x_plot(end) + max(data_x(1,:))/40];
y_plot = (min(data_x(2,:)) : max(data_x(2,:))/40 :max(data_x(2,:)));
y_plot = 1.05*[y_plot, y_plot(end) + max(data_x(2,:))/40];
%
nx_plot = size(x_plot,2);
ny_plot = size(y_plot,2);
for k=1:3
    for i=1:nx_plot
        for j=1:ny_plot
            xx(i,j) = x_plot(i);
            yy(i,j) = y_plot(j);
            plot_xy = [x_plot(i); y_plot(j)];
            vec_a_at_x  = zeros(num.dim,1);
            cons_c_at_x = 0;
            for l=1:num.data
                x_l = data_x(:,l);
                kappa_l = exp(-param.RBF * ((plot_xy-x_l(1:2))'*(plot_xy-x_l(1:2))) );
                vec_a_at_x  = vec_a_at_x  + (kappa_l * vec_a{k,l});
                cons_c_at_x = cons_c_at_x + (kappa_l * cons_c{k,l});
            end
            zz{k}(i,j) = (cons_c_at_x - (vec_a_at_x(1:2)' * plot_xy)) / vec_a_at_x(3);
            clear plot_xy
        end
    end
end
clear x_plot y_plot nx_plot ny_plot
% <---
% Prepare figures: plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
% --->
plot3(data_x(1,:), data_x(2,:), data_x(3,:), 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3);
hold on;
grid on;
xlabel('$x_{1}$', 'Interpreter', 'latex');
ylabel('$x_{2}$', 'Interpreter', 'latex');
zlabel('$y$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
view(-20,25);

figure;
colormap winter;
mesh(xx,yy,zz{1});
hold on;
colormap autumn;
mesh(xx,yy,zz{2});
plot3(data_x(1,:), data_x(2,:), data_x(3,:), 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3);
grid on;
xlabel('$x_{1}$', 'Interpreter', 'latex');
ylabel('$x_{2}$', 'Interpreter', 'latex');
zlabel('$y$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
view(-20,25);


figure;
colormap winter;
mesh(xx,yy,zz{1});
hold on;
colormap autumn;
mesh(xx,yy,zz{2});
plot3(data_x(1,:), data_x(2,:), data_x(3,:), 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3);
grid on;
xlabel('$x_{1}$', 'Interpreter', 'latex');
ylabel('$x_{2}$', 'Interpreter', 'latex');
zlabel('$y$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
view(-5,60);
zlim([-10,10]);
% <---
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensionality detection
% --->
square_error = zeros(1,3);
error_analysis_failed = cell(1,3);

options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,...
    'Display','off');
for k=1:3
    param.error_cur = k;
    for i=1:num.data
        data_cur = data_x(:,i);
        [~,error_cur,exitflag] = fmincon(@error_obj, data_cur, [],[], [],[], [],[],...
            @error_cstr, options);
        if exitflag > 0
            square_error(k) = square_error(k) + error_cur;
        else
            error_analysis_failed{k} = [error_analysis_failed{k}, i];
        end
    end
end

fprintf(' ============================================= \n');
fprintf('   dimensionality analysis \n');
fprintf('      n-1: analytical_error = %3.5d \n',...
    vec_sol{1}' * (matH') * matH * vec_sol{1});
for k=1:3
    fprintf('      n-%g: sum_of_sq_error = %3.5d (failed: %g times) \n',...
        k, square_error(k), length(error_analysis_failed{k}) );
end
fprintf(' ============================================= \n');
% <---
% Dimensionality detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off;

