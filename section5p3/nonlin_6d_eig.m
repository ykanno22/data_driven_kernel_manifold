%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Copyright (C) Yoshihiro Kanno, 2020. All rights reserved %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
% --->
global param num vec_a cons_c data_x data_cur
% <---
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.RBF   =  25.0 * 10^(0);
param.regularize = 1.0 * 10^(-4);
%
load('nonlin_iso_data_set.mat');
num.data = size(list_eps, 2);
num.dim   = 6;
num.dim_x = num.dim / 2;
num.variable = num.data * (num.dim + 1);
%
data_x(1:3,:) = list_eps;
data_x(4:6,:) = list_sigma;
clear list_eps list_sigma
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare matrices & vectors
% - -->
matKappa = zeros(num.data,num.data);
vec_h = cell(1,num.data);
d_x = num.dim_x;
for i=1:num.data
    x_i = data_x(:,i);
    vec_h{i} = [];
    for j=1:num.data
        x_j = data_x(:,j);
        kappa_ij = exp(-param.RBF * ((x_i(1:d_x)-x_j(1:d_x))'*(x_i(1:d_x)-x_j(1:d_x))) );
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
        kappa_ij = exp(-1.0*param.RBF * ((x_i(1:d_x)-x_j(1:d_x))'*(x_i(1:d_x)-x_j(1:d_x))) );
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
vec_sol = cell(1, num.dim);
for k=1:num.dim
    vec_sol{k} = eig_vec(:,k);
end
%
vec_a  = cell(num.dim, num.data);
cons_c = cell(num.dim, num.data);
for k=1:num.dim
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
% Dimensionality analysis
% --->
square_error = zeros(1,num.dim);
error_analysis_failed = cell(1,num.dim);

options = optimoptions('fmincon', 'SpecifyObjectiveGradient',true,...
    'Display','off');
for k=1:num.dim
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
for k=1:num.dim
    fprintf('      n-%g: sum_of_sq_error = %3.5d (failed: %g times) \n',...
        k, square_error(k), length(error_analysis_failed{k}) );
end
fprintf(' ============================================= \n');
% <---
% Dimensionality analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num.plane = 3;
np = num.plane;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('nonlin_eig_6d.mat',...
    'param', 'num', 'vec_a', 'cons_c', 'data_x', 'param_scale',...
    'mean_Young', 'mean_Poisson',...
    'list_norm_d_eps', 'list_norm_d_sig');

