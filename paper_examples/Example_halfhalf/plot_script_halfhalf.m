% Script for generating Figure 2(a) in the article 

clear all
rng('default')

addpath(genpath(fullfile('../../test_funs/halfhalf/')));

n = 8;
problem_data.n = n;
problem_data.f = @halfhalf;
problem_data.subgrad_f = @grad_halfhalf;
problem_data.subhess_f = @hess_halfhalf;
problem_data.x0 = 20.08*ones(n,1);

%% Set parameters for SOGS

j_max = 5;
eps0 = 10;
kappa_eps = 0.1;
eps_fun = @(j) eps0*kappa_eps.^(j-1);
sogs_options.eps_arr = eps_fun(1:j_max);
sogs_options.tau_arr = 10^-5*ones(1,j_max);

sogs_options.c = 0.5;
sogs_options.sp_solver = 'IPOPT';
sogs_options.sp_solver_options.tol = 10^-10;
% sogs_options.sp_solver = 'fmincon';
% sogs_options.sp_solver_options = optimoptions('fmincon');
% sogs_options.sp_solver_options.Display = 'none';
% sogs_options.sp_solver_options.Algorithm = 'active-set';
% sogs_options.sp_solver_options.OptimalityTolerance = 10^-10;
% sogs_options.sp_solver_options.FunctionTolerance = 10^-10;
% sogs_options.sp_solver_options.StepTolerance = 10^-10;
sogs_options.max_iter = 10000;
sogs_options.init_N_sample = 1;
sogs_options.memory_max_size = 100;
sogs_options.disp_flag = 1;

%% Run the algorithm

addpath(genpath(fullfile('../..')));
[~,~,x_cell,~,numsample_cell,mu_cell] = sogs(problem_data,sogs_options);

%% Process the results
f = problem_data.f;

xj_arr = zeros(n,j_max);
fxj_arr = zeros(1,j_max);
xji_arr = [];
ji_arr = NaN(1,j_max);
numsamplej_arr = zeros(1,j_max);
numsampleji_arr = [];
for j = 1:j_max
    xj_arr(:,j) = x_cell{j}(:,end);
    fxj_arr(j) = f(x_cell{j}(:,end));
    xji_arr = [xji_arr,x_cell{j}];
    ji_arr(j) = size(xji_arr,2);
    numsamplej_arr(j) = sum(numsample_cell{j});
    numsampleji_arr = [numsampleji_arr,numsample_cell{j}];

    % Append the final point to include all evaluations. (Otherwise,
    % evaluations for the stopping criterion in the final iterate would be
    % ignored.)
    if(j == j_max)
        xji_arr = [xji_arr,x_cell{j}(:,end)];
        numsampleji_arr = [numsampleji_arr,0];
    end
end
l_max = size(xji_arr,2);

fxji_arr = zeros(1,l_max);
for l = 1:l_max
    fxji_arr(l) = f(xji_arr(:,l));
end

%% Visualization

% Load result of VU_bundle
load("julia-results-03-01-25_10-25-50")

x_min = zeros(n,1);
f_min = f(x_min);

lw = 1.5;
ms = 14;

subgrad_eval_arr = tril(ones(l_max) - eye(l_max))*numsampleji_arr';

figure
h1 = plot(subgrad_eval_arr,log10(vecnorm(fxji_arr - f_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on
h2 = plot(subgrad_eval_arr(ji_arr),log10(vecnorm(fxji_arr(ji_arr) - f_min,2,1)),'r.','LineWidth',lw,'MarkerSize',ms+0.5);
plot(subgrad_eval_arr([1,ji_arr]),log10(vecnorm(fxji_arr([1,ji_arr]) - f_min,2,1)),'r--','LineWidth',lw,'MarkerSize',ms+0.5);
h3 = plot(subgrad_arr_vu,log10(f_arr_vu - f_min),'b.-','LineWidth',lw,'MarkerSize',ms);
ylabel('Distance to optimal value','Interpreter','latex');
grid on
xlabel('Subgradient evaluations','Interpreter','latex');

box on
legend([h1,h3],{'\verb|SOGS|','\verb|VUbundle|'},'Interpreter','latex','Location','northeast','FontSize',20)

ylim([-11,4])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

return
%%

export_fig 'plot' '-jpg' '-r500' '-transparent'
