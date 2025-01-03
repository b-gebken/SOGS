% This script contains a simple example showing how to use SOGS to solve
% the Half-and-half problem from [1]. 

% [1] Mifflin, Sagastizábal, "A science fiction story in nonsmooth
% optimization originating at IIASA" (2012)

%% Define the problem

addpath('test_funs/halfhalf');
n = 8;
problem_data.n = n;
problem_data.f = @halfhalf;
problem_data.subgrad_f = @grad_halfhalf;
problem_data.subhess_f = @hess_halfhalf;
problem_data.x0 = 20.08*ones(n,1);

%% Set the parameters for SOGS

% Generate tolerance arrays eps_j and tau_j
j_max = 5;
eps0 = 10;
kappa_eps = 0.1;
eps_fun = @(j) eps0*kappa_eps.^(j-1);
sogs_options.eps_arr = eps_fun(1:j_max);
sogs_options.tau_arr = 10^-5*ones(1,j_max);

% Choose solver and set the options (IPOPT is significantly faster than
% fmincon in our experience) 
sogs_options.sp_solver = 'IPOPT';
sogs_options.sp_solver_options.tol = 10^-10;

% sogs_options.sp_solver = 'fmincon';
% sogs_options.sp_solver_options = optimoptions('fmincon');
% sogs_options.sp_solver_options.Display = 'none';
% sogs_options.sp_solver_options.Algorithm = 'active-set';
% sogs_options.sp_solver_options.OptimalityTolerance = 10^-10;
% sogs_options.sp_solver_options.FunctionTolerance = 10^-10;
% sogs_options.sp_solver_options.StepTolerance = 10^-10;

% Controls the quality of the jet approximation
sogs_options.c = 0.5;

% The number of randomly sampled points (1 means no random sampling) 
sogs_options.init_N_sample = 1;

% The maximum number of jet elements that are carried over between
% iterations
sogs_options.memory_max_size = 100;

% Maximum number of (inner) iterations
sogs_options.max_iter = 1000; 

% Flag for controlling the amount of output that is displayed during the
% algorithm 
sogs_options.disp_flag = 1; 

%% Run the algorithm

[x_opt,f_opt,x_cell,eval_counter,numsample_cell,mu_cell] = sogs(problem_data,sogs_options);

%% Visualization

x_min = zeros(n,1);
f_min = problem_data.f(x_min);

xl_arr = cell2mat(x_cell);
l_max = size(xl_arr,2);
fxl_arr = zeros(1,l_max);
for i = 1:l_max
    fxl_arr(i) = problem_data.f(xl_arr(:,i));
end

lw = 1.5;
ms = 12;

numsamplel_arr = cell2mat(numsample_cell); 
subgrad_eval_arr = tril(ones(l_max) - eye(l_max))*numsamplel_arr';

figure
subplot(1,2,1)
h1 = plot(1:l_max,log10(vecnorm(fxl_arr - f_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
ylabel('Distance to optimal value','Interpreter','latex');
grid on
xlabel('Iterations','Interpreter','latex');

legend(h1,{'\verb|SOGS|'},'Interpreter','latex','Location','northeast','FontSize',20)
xlim([0,l_max+1]);
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

subplot(1,2,2);
h2 = plot(subgrad_eval_arr,log10(vecnorm(fxl_arr - f_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
ylabel('Distance to optimal value','Interpreter','latex');
grid on
xlabel('Subgradient evaluations','Interpreter','latex');

legend(h2,{'\verb|SOGS|'},'Interpreter','latex','Location','northeast','FontSize',20)
xlim([0,subgrad_eval_arr(end)+1])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)