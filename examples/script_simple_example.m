%% This script shows a simple application of our method.

% Define the problem (in this case "Crescent" from Kiwiel (1985))
addpath(fullfile(cd,'..','test_funs\H2004lib\Chained_Crescent_I'));
problem_config.n = 2;
problem_config.x0 = [-1.5;2];
problem_config.f = @(x) crescent_I(x);
problem_config.subgrad_f = @(x) grad_crescent_I(x);
problem_config.subhess_f = @(x) hess_crescent_I(x);

% Set the parameters for the method (see "second_order_gradient_sampling.m" for an explanation)
SOGS_options.eps_init = 1*10^1;
SOGS_options.kappa_eps = 10^-1;
SOGS_options.eps_end = 1*10^-5;

SOGS_options.tau_init = 1*10^-5;
SOGS_options.tau_end = 1*10^-5;

SOGS_options.reuse_flag = 1;
SOGS_options.init_N_sample = 1;
SOGS_options.c = 0.5;
SOGS_options.max_iter = 1000;
SOGS_options.disp_flag = 2;

%% Run the method
addpath(fullfile(cd,'..')); % Path of the second-order gradient sampling method

[x_arr_SOGS,fval,numsample_arr_SOGS,eps_ind_SOGS,lambda_arr] = second_order_gradient_sampling(problem_config,SOGS_options);
sol = x_arr_SOGS(:,end);

%% Plot the result
[X1,X2] = ndgrid(-2:0.05:2);
fX = zeros(size(X1));
for i = 1:size(X1,1)
    for j = 1:size(X1,2)
        fX(i,j) = problem_config.f([X1(i,j);X2(i,j)]);
    end
end

figure
surf(X1,X2,fX);
hold on
plot3(sol(1),sol(2),fval,'r.','MarkerSize',30);
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$f$','Interpreter','latex')

