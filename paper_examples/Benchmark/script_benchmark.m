% A script that applies all six solvers to the 20 test problems. For a fair
% comparison, the number of function, subgradient and Hessian evaluations
% are counted in the same way for each solver (via global variables).
% Furthermore, there are sections for saving and analyzing the results.
% (These can be executed via a Ctrl+Enter.)

clear all
rng('default')

%% Define the test problems
n = 50;
problem_inds = 1:20;
num_problems = numel(problem_inds);

problem_arr(num_problems) = struct();
addpath(genpath(fullfile('../../test_funs/H2004lib/')));
addpath(genpath(fullfile('../../test_funs/TEST29/')));
func_name_cell = {'maxq','mxhilb','chained_LQ','chained_CB3_I','chained_CB3_II','active_faces','brown_2','chained_mifflin_2','crescent_I','crescent_II'};
TEST29_ind_shift = [2,5,6,11,13,17,19,20,22,24];

for k = 1:num_problems
    prob_i = problem_inds(k);
    problem_arr(k).n = n;

    if(prob_i <= 10)
        % HMM2004
        if(prob_i == 7)
            eval(['init_',func_name_cell{prob_i}]);
            problem_arr(k).f = @(in) counter_wrapper(brown_2(in),1);
            problem_arr(k).subgrad_f = @(in) counter_wrapper(grad_brown_2(eval_pts_wrapper(in)),2);
            problem_arr(k).subhess_f = @(in) counter_wrapper(hess_brown_2(in),3);
        else
            problem_arr(k).f = @(in) counter_wrapper(feval(func_name_cell{prob_i},in),1);
            problem_arr(k).subgrad_f = @(in) counter_wrapper(feval(['grad_',func_name_cell{prob_i}],eval_pts_wrapper(in)),2);
            problem_arr(k).subhess_f = @(in) counter_wrapper(feval(['hess_',func_name_cell{prob_i}],in),3);
        end

        problem_arr(k).x0 = feval([func_name_cell{prob_i},'_x0'],n);
    else
        % TEST29 from H2004
        problem_arr(k).f = @(in) counter_wrapper(feval(['problem',num2str(TEST29_ind_shift(prob_i - 10))],in),1);
        problem_arr(k).subgrad_f = @(in) counter_wrapper(feval(['grad_problem',num2str(TEST29_ind_shift(prob_i - 10))],eval_pts_wrapper(in)),2);
        problem_arr(k).subhess_f = @(in) counter_wrapper(feval(['hess_problem',num2str(TEST29_ind_shift(prob_i - 10))],in),3);

        problem_arr(k).x0 = feval(['problem',num2str(TEST29_ind_shift(prob_i - 10)),'_x0'],n);
    end

end

%% Set parameters for all solvers

% SOGS parameters
j_max = 5;
eps0 = 10;
kappa_eps = 0.1;
eps_fun = @(j) eps0*kappa_eps.^(j-1);
sogs_options.eps_arr = eps_fun(1:j_max);
sogs_options.tau_arr = 10^-5*ones(1,j_max);
sogs_options.c = 0.5;
sogs_options.sp_solver = 'IPOPT';
sogs_options.sp_solver_options.tol = 10^-8;
sogs_options.max_iter = 10000; % (This is never reached.)
sogs_options.init_N_sample = 1;
sogs_options.memory_max_size = 100;
sogs_options.disp_flag = 0;

% DGS parameters
j_max = 7;
kappa_eps = 0.1;
eps0 = 10;
eps_fun = @(j) eps0*kappa_eps.^(j);
del_fun = @(j) 10^-3 * ones(size(j));
dgs_options.eps_arr = eps_fun(0:j_max-1);
dgs_options.delta_arr = del_fun(0:j_max-1);
dgs_options.rand_sample_N = 0;
dgs_options.memory_size = Inf;
dgs_options.c = 0.5;
dgs_options.ls_flag = 'armijo_normal';
dgs_options.max_iter = sogs_options.max_iter;
dgs_options.disp_flag = 0;

% Gradsamp parameters
gradsamp_options.prtlevel = 0;
gradsamp_options.maxit = sogs_options.max_iter;

% Hanso parameters
hanso_options.prtlevel = 0;
hanso_options.maxit = sogs_options.max_iter;

% SLQPGS parameters
slqpgs_options.iter_max = sogs_options.max_iter;
slqpgs_options.nE = 0;
slqpgs_options.nI = 0;
slqpgs_options.pE = 0;
slqpgs_options.pI = 0;
slqpgs_options.f = 'my_slqpgs_problem';
slqpgs_options.output = 0;
slqpgs_options.log_fields = {'x'};

% LMBM parameters (for default values, see H2004, p. 85 and lmbm.f, line
% 574) 
lmbm_options.maxtime = 300;
lmbm_options.print = 0;
lmbm_options.na = 100; % See Section 6.2.2 in H2004
lmbm_options.mcu = 50; % See p. 90 in H2004
lmbm_options.mc = 7; % See H2004
lmbm_options.rpar = [10^-8, 10^4, -10^60, 10^-6, 10^-6, 0.5, 10^-4, 1.5];
lmbm_options.ipar = [2, sogs_options.max_iter, 20000, 10, -1, 0, 0];

%% Solve all problems with all solvers

% Preallocate structs containing the output of the solvers
sogs_output = repmat(struct([]),[1,num_problems]);
gradsamp_output = repmat(struct([]),[1,num_problems]);
dgs_output = repmat(struct([]),[1,num_problems]);
hanso_output = repmat(struct([]),[1,num_problems]);
slqpgs_output = repmat(struct([]),[1,num_problems]);
lmbm_output = repmat(struct([]),[1,num_problems]);

% Create global variables
global f_counter
global subgrad_counter
global subhess_counter

global eval_pts_global % For recording subgrad. eval. points

global function_data % For LMBM

totalT = tic;
for i = 1:num_problems
    disp(['----- Problem ',num2str(problem_inds(i)),' -----'])

    % SOGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];
    addpath('../..');
    disp('SOGS...')
    tic
    [~,f_opt,x_cell,~,numsample_cell,mu_cell] = sogs(problem_arr(i),sogs_options);
    sogs_output(i).runtime = toc;
    disp('    ...done!')
    rmpath('../..');
    sogs_output(i).eval_pts = eval_pts_global;
    sogs_output(i).x_cell = x_cell;
    sogs_output(i).x_arr = [x_cell{:}];
    sogs_output(i).num_iter = size(sogs_output(i).x_arr,2)-1;
    sogs_output(i).numsample_cell = numsample_cell;
    sogs_output(i).mu_cell = mu_cell;
    sogs_output(i).f_eval = f_counter;
    sogs_output(i).subgrad_eval = subgrad_counter;
    sogs_output(i).subhess_eval = subhess_counter;
    sogs_output(i).final_f = f_opt;
    
    % Gradsamp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pars.m = problem_arr(i).n;
    pars.fgname = 'my_fun';
    pars.my_f = problem_arr(i).f;
    pars.my_subgrad_f = problem_arr(i).subgrad_f;
    gradsamp_options.x0 = problem_arr(i).x0;
    
    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];
    addpath('gradsamp');
    disp('Gradsamp...')
    tic
    [~,f_opt,optcert] = gradsamp(pars,gradsamp_options);
    gradsamp_output(i).runtime = toc;
    disp('    ...done!')
    rmpath('gradsamp');
    gradsamp_output(i).eval_pts = eval_pts_global;
    gradsamp_output(i).f_eval = f_counter;
    gradsamp_output(i).subgrad_eval = subgrad_counter;
    gradsamp_output(i).subhess_eval = subhess_counter;
    gradsamp_output(i).final_f = f_opt;
    gradsamp_output(i).optcert = optcert;
    clear pars
    
    % DGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];
    addpath('DGS-main');
    disp('DGS...')
    tic
    [~,f_opt,x_cell,~] = eps_descent_method(problem_arr(i),dgs_options);
    dgs_output(i).runtime = toc;
    disp('    ...done!')
    rmpath('DGS-main');
    dgs_output(i).eval_pts = eval_pts_global;
    dgs_output(i).x_cell = x_cell;
    dgs_output(i).x_arr = [x_cell{:}];
    dgs_output(i).num_iter = size(dgs_output(i).x_arr,2)-1;
    dgs_output(i).f_eval = f_counter;
    dgs_output(i).subgrad_eval = subgrad_counter;
    dgs_output(i).subhess_eval = subhess_counter;
    dgs_output(i).final_f = f_opt;

    % Hanso %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pars.nvar = problem_arr(i).n;
    pars.fgname = 'my_fun';
    pars.my_f = problem_arr(i).f;
    pars.my_subgrad_f = problem_arr(i).subgrad_f;
    hanso_options.x0 = problem_arr(i).x0;

    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];
    addpath('hanso3_0');
    disp('Hanso...')
    tic
    [~,f_opt,loc] = hanso(pars, hanso_options);
    hanso_output(i).runtime = toc;
    disp('    ...done!')
    rmpath('hanso3_0');
    hanso_output(i).eval_pts = eval_pts_global;
    hanso_output(i).f_eval = f_counter;
    hanso_output(i).subgrad_eval = subgrad_counter;
    hanso_output(i).subhess_eval = subhess_counter;
    hanso_output(i).final_f = f_opt;
    hanso_output(i).loc = loc;
    clear pars

    % SLQPGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slqpgs_options.nV = problem_arr(i).n;
    slqpgs_options.pO = 2*problem_arr(i).n;
    slqpgs_options.x = problem_arr(i).x0;
    slqpgs_options.d.f = problem_arr(i).f;
    slqpgs_options.d.subgrad_f = problem_arr(i).subgrad_f;
    
    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];
    addpath('SLQPGS/src/');
    slqpgs_prob = SLQPGS(slqpgs_options);
    disp('SLQPGS...')
    tic
    slqpgs_prob.optimize;
    slqpgs_output(i).runtime = toc;
    disp('    ...done!')
    [x_opt,slqpgs_log] = slqpgs_prob.getSolution();
    clear slqpgs_prob
    rmpath('SLQPGS/src/');
    slqpgs_output(i).eval_pts = eval_pts_global;
    tmp_cell = cell(1,size(slqpgs_log,2));
    [tmp_cell{:}] = deal(slqpgs_log.x);
    slqpgs_output(i).x_arr = cell2mat(tmp_cell);
    slqpgs_output(i).num_iter = size(slqpgs_output(i).x_arr,2)-1;
    slqpgs_output(i).f_eval = f_counter;
    slqpgs_output(i).subgrad_eval = subgrad_counter;
    slqpgs_output(i).subhess_eval = subhess_counter;
    slqpgs_output(i).final_f = problem_arr(i).f(x_opt);
    slqpgs_output(i).log = slqpgs_log;

    % LMBM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_counter = 0; subgrad_counter = 0; subhess_counter = 0;
    eval_pts_global = [];

    addpath('lmbm-mex');
    function_data.f = @(in) problem_arr(i).f(lmbm_eval_wrapper(in));
    function_data.subgrad_f = problem_arr(i).subgrad_f;
    disp('LMBM...')
    tic
    % (1* is added because pointers in C)
    [~,fval,niter,nfeval,term]=lmbm_driver('my_fun', 1*problem_arr(i).x0, 1*problem_arr(i).n, lmbm_options.print,...
        lmbm_options.maxtime, lmbm_options.na, lmbm_options.mcu, lmbm_options.mc, lmbm_options.rpar, lmbm_options.ipar);
    lmbm_output(i).runtime = toc;
    disp('    ...done!')
    rmpath('lmbm-mex')
    lmbm_output(i).eval_pts = eval_pts_global;
    lmbm_output(i).f_eval = nfeval;
    lmbm_output(i).subgrad_eval = nfeval;
    lmbm_output(i).subhess_eval = 0;
    lmbm_output(i).final_f = fval;
    lmbm_output(i).niter = niter;
end
bench_time = toc(totalT);
disp(['Benchmark completed in ',num2str(bench_time/60),' minutes.'])

return
%% Save results

save(['results_',datestr(now,'dd-mm-yy_HH-MM-SS')],'problem_arr',...
    'sogs_output','gradsamp_output','dgs_output','hanso_output','slqpgs_output','lmbm_output',...
    'sogs_options','gradsamp_options','dgs_options','hanso_options','slqpgs_options','lmbm_options');

%% Present the results as tables

load('f_min_arr.mat');

for i = 1:numel(problem_arr)
    disp(['----- Problem ',num2str(i),' -----'])
    perf_table(f_min_arr(i),sogs_output(i),gradsamp_output(i),dgs_output(i),hanso_output(i),slqpgs_output(i),lmbm_output(i));
end

%% Plots for specific problems (chosen via k below)

addpath(genpath(fullfile('../../test_funs/H2004lib/')));
addpath(genpath(fullfile('../../test_funs/TEST29/')));

k = 1;
n = problem_arr(k).n;
conv_thresh = 10^-4;

load('f_min_arr')
f_min = f_min_arr(k);

% Compute function value arrays
f_arr_sogs = zeros(1,size(sogs_output(k).eval_pts,2));
for i = 1:size(sogs_output(k).eval_pts,2)
    f_arr_sogs(i) = problem_arr(k).f(sogs_output(k).eval_pts(:,i)); 
end

f_arr_gradsamp = zeros(1,size(gradsamp_output(k).eval_pts,2));
for i = 1:size(gradsamp_output(k).eval_pts,2)
    f_arr_gradsamp(i) = problem_arr(k).f(gradsamp_output(k).eval_pts(:,i)); 
end

f_arr_dgs = zeros(1,size(dgs_output(k).eval_pts,2));
for i = 1:size(dgs_output(k).eval_pts,2)
    f_arr_dgs(i) = problem_arr(k).f(dgs_output(k).eval_pts(:,i)); 
end

f_arr_hanso = zeros(1,size(hanso_output(k).eval_pts,2));
for i = 1:size(hanso_output(k).eval_pts,2)
    f_arr_hanso(i) = problem_arr(k).f(hanso_output(k).eval_pts(:,i)); 
end

f_arr_slqpgs = zeros(1,size(slqpgs_output(k).eval_pts,2));
for i = 1:size(slqpgs_output(k).eval_pts,2)
    f_arr_slqpgs(i) = problem_arr(k).f(slqpgs_output(k).eval_pts(:,i)); 
end

f_arr_lmbm = zeros(1,size(lmbm_output(k).eval_pts,2));
for i = 1:size(lmbm_output(k).eval_pts,2)
    f_arr_lmbm(i) = problem_arr(k).f(lmbm_output(k).eval_pts(:,i)); 
end

disp(['Problem ',num2str(k)]);
disp('Cutoffs');
disp('    SOGS')
disp(['        #subgrads: ',num2str(find(f_arr_sogs - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_sogs(find(f_arr_sogs - f_min < conv_thresh,1,'first')) - f_min)])
disp('    Gradsamp')
disp(['        #subgrads: ',num2str(find(f_arr_gradsamp - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_gradsamp(find(f_arr_gradsamp - f_min < conv_thresh,1,'first')) - f_min)])
disp('    DGS')
disp(['        #subgrads: ',num2str(find(f_arr_dgs - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_dgs(find(f_arr_dgs - f_min < conv_thresh,1,'first')) - f_min)])
disp('    Hanso')
disp(['        #subgrads: ',num2str(find(f_arr_hanso - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_hanso(find(f_arr_hanso - f_min < conv_thresh,1,'first')) - f_min)])
disp('    SLQPGS')
disp(['        #subgrads: ',num2str(find(f_arr_slqpgs - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_slqpgs(find(f_arr_slqpgs - f_min < conv_thresh,1,'first')) - f_min)])
disp('    LMBM')
disp(['        #subgrads: ',num2str(find(f_arr_lmbm - f_min < conv_thresh,1,'first'))])
disp(['        f:         ',num2str(f_arr_lmbm(find(f_arr_lmbm - f_min < conv_thresh,1,'first')) - f_min)])

figure
h1 = plot(1:size(f_arr_sogs,2),log10(f_arr_sogs - f_min),'.-');
hold on
h2 = plot(1:size(f_arr_gradsamp,2),log10(f_arr_gradsamp - f_min),'.-');
h3 = plot(1:size(f_arr_dgs,2),log10(f_arr_dgs - f_min),'.-');
h4 = plot(1:size(f_arr_hanso,2),log10(f_arr_hanso - f_min),'.-');
h5 = plot(1:size(f_arr_slqpgs,2),log10(f_arr_slqpgs - f_min),'.-');
h6 = plot(1:size(f_arr_lmbm,2),log10(f_arr_lmbm - f_min),'.-');
legend([h1,h2,h3,h4,h5,h6],{'SOGS','gradsamp','DGS','Hanso','SLQPGS','LMBM'});
yline(log10(conv_thresh),'k--')
title(['Problem ',num2str(k)])
grid on

%% An array containing the best function values be encountered 

n = 50;
f_min_arr = [0; % 1
    0; % 2
    -(n-1)*2^(1/2); % 3
    2*(n-1); % 4
    2*(n-1); % 5
    0; % 6
    0; % 7
    -34.795181336028939; % 8
    0; % 9
    0; % 10
    0; % 11
    0; % 12
    0; % 13
    58.1845604639973935983977516051; % 14
    27.227867622575232; % 15
    0; % 16
    0; % 17
    0; % 18
    0; % 19
    0]; % 20

save('f_min_arr','f_min_arr');