%% This script reproduces the results in Table 1.

% Define problems
problem_arr = struct('n', repmat({[]}, 1, 20),'x0',[],'f',[],'subgrad_f',[],'subhess_f',[]);

% HMM2004
addpath(genpath(fullfile(cd,'..','test_funs/H2004lib/')));
func_name_cell = {'maxq','mxhilb','chained_LQ','chained_CB3_I','chained_CB3_II','active_faces','brown_2','mifflin_2','crescent_I','crescent_II'};
for prob_no = 1:10
    problem_arr(prob_no).n = 50;
    if(prob_no == 7 || prob_no == 8)
        n = problem_arr(prob_no).n;
        eval(['init_',func_name_cell{prob_no}])
        eval(['f = ',func_name_cell{prob_no},'; subgrad_f = grad_',func_name_cell{prob_no},'; subhess_f = hess_',func_name_cell{prob_no},';']);
        problem_arr(prob_no).f = f;
        problem_arr(prob_no).subgrad_f = subgrad_f;
        problem_arr(prob_no).subhess_f = subhess_f;
        clear f subgrad_f subhess_f
    else
        problem_arr(prob_no).f = @(in) feval(func_name_cell{prob_no},in);
        problem_arr(prob_no).subgrad_f = @(in) feval(['grad_',func_name_cell{prob_no}],in);
        problem_arr(prob_no).subhess_f = @(in) feval(['hess_',func_name_cell{prob_no}],in);
    end

    problem_arr(prob_no).x0 = feval([func_name_cell{prob_no},'_x0'],problem_arr(prob_no).n);
end

% TEST29 from H2004
addpath(genpath(fullfile(cd,'..','test_funs/TEST29/')));
prob_no_shift = [2,5,6,11,13,17,19,20,22,24];

for prob_no = 1:10
    problem_arr(prob_no+10).n = 50;

    problem_arr(prob_no+10).f = @(in) feval(['problem',num2str(prob_no_shift(prob_no))],in);
    problem_arr(prob_no+10).subgrad_f = @(in) feval(['grad_problem',num2str(prob_no_shift(prob_no))],in);
    problem_arr(prob_no+10).subhess_f = @(in) feval(['hess_problem',num2str(prob_no_shift(prob_no))],in);

    problem_arr(prob_no+10).x0 = feval(['problem',num2str(prob_no_shift(prob_no)),'_x0'],problem_arr(prob_no+10).n);
end

num_problems = numel(problem_arr);

%% Set parameters for algorithms

% SOGS config
SOGS_options.eps_init = 1*10^1;
SOGS_options.kappa_eps = 10^-1;
SOGS_options.eps_end = 1*10^-5;

SOGS_options.tau_init = 1*10^-5;
SOGS_options.tau_end = 1*10^-5;

SOGS_options.reuse_flag = 1;
SOGS_options.init_N_sample = 1;
SOGS_options.c = 0.5;
SOGS_options.max_iter = 1000;
SOGS_options.disp_flag = 1;

%% Run the benchmark

addpath(fullfile(cd,'..')); % Path of the second-order gradient sampling method

% SOGS output
SOGS_output = struct('x_arr', repmat({[]}, 1, num_problems),'numsample_arr',[],'eps_ind',[],'lambda_arr',[],'runtime',[],'final_f',[]);

totalT = tic;
for i = 1:num_problems
    disp(['----- Problem ',num2str(i),' -----'])
    
    disp('Second-order gradient sampling...')
    tic
    [x_arr_SOGS,fval,numsample_arr_SOGS,eps_ind_SOGS,lambda_arr] = second_order_gradient_sampling(problem_arr(i),SOGS_options);
    SOGS_output(i).runtime = toc;
    disp('    ...done!')
    SOGS_output(i).x_arr = x_arr_SOGS;
    SOGS_output(i).numsample_arr = numsample_arr_SOGS;
    SOGS_output(i).eps_ind = eps_ind_SOGS;
    SOGS_output(i).lambda_arr = lambda_arr;
    SOGS_output(i).final_f = fval;
    
end
bench_time = toc(totalT);
disp(['Benchmark completed in ',num2str(bench_time/60),' minutes.'])

%% Display the result

subgrad_eval = cellfun(@(in) sum(in),{SOGS_output.numsample_arr});
results_table = table((1:num_problems)',subgrad_eval',[SOGS_output.final_f]','VariableNames',{'Problem','Evaluations','Final_value'});
disp(results_table)
