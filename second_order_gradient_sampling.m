function [x_arr,fval,numsample_arr,eps_ind,lambda_arr] = second_order_gradient_sampling(problem_config,algo_options)
%   Second-order gradient sampling method for nonsmooth optimization
%       An implementation of Algorithm 3 from [1] for the minimization of
%       nonsmooth, nonconvex optimization problems based on sampling 
%       second-order information of the objective function.
%       
%   Input:
%       problem_config is a struct with the following fields:
%           n: Number of variables.
%           x0: Initial point (as a column vector).
%           f: Function handle that returns the objective function at a
%               given point.
%           subgrad_f: Function handle that returns an arbitrary element
%               of the subdifferential at a given point.
%           subhess_f: Function handle that returns a matrix for the
%               second-order eps-jet at a given point. (Typically the
%               Hessian matrix in that point.)
%
%       algo_options is a struct with the followings fields:
%           eps_init: Initial value of epsilon.
%           kappa_eps: Reduction factor for epsilon.
%           eps_end: Final value of epsilon (for stopping).
%           tau_init: Initial value of tau.
%           tau_end: Final value of epsilon. (The reduction factor
%               kappa_tau in [1] is computed from the required reduction 
%               steps for epsilon.)        
%           reuse_flag: A flag for the initialization of W in step 2:
%               0: The initial W only contains an element of the 0-jet at
%               the current point x^i (as in step 1 of Algo. 2).
%               1: Already sampled elements of the eps-jet from the
%               previous iteration that also lie in the current eps-jet are
%               reused.
%           init_N_sample: Number of randomly sampled points (in addition 
%               to potentially reused points) for the initialization of W
%               in step 2. If init_N_sample = 1, then only an element at 
%               the current point x^i is used (as in step 1 of Algo. 2).
%           c: Approximation parameter for the eps-jet (in (0,1)).
%           max_iter: Maximum number of iterations.
%           disp_flag: A flag for the output in the command window during
%               the algorithm:
%               0: No output (except for possible warnings)
%               1: Final objective value and number of subgradient 
%               evaluations are displayed after the run.
%               2: Details on every iteration of the algorithm. This
%               includes the norm of the element with the smallest norm in
%               the convex hull of all evaluated subgradients (computed via 
%               quadprog), which can be used as a measure for criticality.
%               (Note that this measure is only displayed and never used in
%               the algorithm.)
%               
%   Output:
%       x_arr: Sequence of generated iterates.
%       fval: Value of the objective in the final iterate.
%       numsample_arr: Number of sampled elements of the eps-jet in each 
%           iteration.
%       eps_ind: Indices of the iterations in which eps was reduced.
%       lambda_arr: Lagrange multiplier of the eps-ball constraint in the
%           subproblem in each iteration.
%
%   [1] Gebken, "Using second-order information in gradient sampling
%   methods for nonsmooth optimization" (2022)

    % Read inputs
    n = problem_config.n;
    x0 = problem_config.x0;
    f = problem_config.f;
    subgrad_f = problem_config.subgrad_f;
    subhess_f = problem_config.subhess_f;
    
    eps_init = algo_options.eps_init;
    kappa_eps = algo_options.kappa_eps;
    eps_end = algo_options.eps_end;
    tau_init = algo_options.tau_init;
    tau_end = algo_options.tau_end;
    c = algo_options.c;
    max_iter = algo_options.max_iter;
    init_N_sample = algo_options.init_N_sample;
    reuse_flag = algo_options.reuse_flag;
    disp_flag = algo_options.disp_flag;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Initialization %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set the options for solving the subproblem via fmincon
    optns = optimoptions('fmincon','Display','off','Algorithm','interior-point');
    optns.MaxIterations = 10^5;
    optns.MaxFunctionEvaluations = 10^5;
    
    x_arr = zeros(n,max_iter+1);
    x_arr(:,1) = x0;
    f_x = f(x0);
    
    numsample_arr = zeros(1,max_iter);
    lambda_arr = zeros(1,max_iter);
    
    kappa_tau = (tau_end/tau_init)^(-1/(log(eps_end/eps_init)/log(1/kappa_eps)));
    num_decr_steps = double(ceil(sym(log(eps_end/eps_init)/log(kappa_eps)))); % Symbolic computation to avoid numerical errors while rounding
    eps = eps_init;
    tau = tau_init;
    eps_ind = zeros(1,num_decr_steps);
    
    sample_pts = [];
    im_f = [];
    subdiff_f = [];
    subhess_f_cell = {};
    
    i = 1;
    j_eps = 0;
    while(j_eps <= num_decr_steps && i <= max_iter)
        if(disp_flag > 1)
            disp(['----- Iteration ',num2str(i),' -------------------']);
            disp(['    eps = ',num2str(eps),', tau = ',num2str(tau)]);
        end
        
        %%%%%%%%%%%%%%%%%%
        %%%%% Step 2 %%%%%
        %%%%%%%%%%%%%%%%%%
        
        % Sample new elements of the eps-jet
        if(disp_flag > 1)
            disp(['    Initially sampling ',num2str(init_N_sample),' point(s).'])
        end
        new_sample_pts = x_arr(:,i) + eps*sample_hypersphere(n,init_N_sample);
        new_im_f = zeros(1,init_N_sample);
        new_subdiff_f = zeros(n,init_N_sample);
        new_subhess_f_cell = cell(1,init_N_sample);
        new_im_f(1) = f_x; % First point is the center of the eps-ball
        for j = 1:init_N_sample
            if(j > 1)
                new_im_f(j) = f(new_sample_pts(:,j));
            end
            new_subdiff_f(:,j) = subgrad_f(new_sample_pts(:,j));
            new_subhess_f_cell{j} = subhess_f(new_sample_pts(:,j));
        end
        
        % Reuse old elements
        if(reuse_flag && size(sample_pts,2) > 0)
            inds = vecnorm(sample_pts - x_arr(:,i),2,1) <= eps;
            sample_pts = sample_pts(:,inds);
            im_f = im_f(inds);
            subdiff_f = subdiff_f(:,inds);
            subhess_f_cell = subhess_f_cell(inds);
            N_reused = sum(inds);
            
            sample_pts = [sample_pts,new_sample_pts];
            im_f = [im_f,new_im_f];
            subdiff_f = [subdiff_f,new_subdiff_f];
            subhess_f_cell = [subhess_f_cell(:)',new_subhess_f_cell(:)'];
            if(disp_flag > 1)
                disp(['    Reusing ',num2str(N_reused),' sample points.'])
            end
        else
            N_reused = 0;
            sample_pts = new_sample_pts;
            im_f = new_im_f;
            subdiff_f = new_subdiff_f;
            subhess_f_cell = new_subhess_f_cell;
        end
        N_sample = size(sample_pts,2);
        
        %%%%%%%%%%%%%%%%%%
        %%%%% Step 3 %%%%%
        %%%%%%%%%%%%%%%%%%
        
        [z,theta,lambda] = solve_subproblem(sample_pts,im_f,subdiff_f,subhess_f_cell,x_arr(:,i),eps,optns);
        f_z = f(z);
        
        %%%%%%%%%%%%%%%%%%
        %%%%% Step 4 %%%%%
        %%%%%%%%%%%%%%%%%%
        
        reduction_flag = false;
        if(theta - f_x > -tau*eps)
            reduction_flag = true;
        end
        
        while(f_z > f_x + c*(theta - f_x) && ~reduction_flag)
           
            %%%%%%%%%%%%%%%%%%
            %%%%% Step 5 %%%%%
            %%%%%%%%%%%%%%%%%%
            
            N_sample = N_sample+1;
            sample_pts = [sample_pts,z];
            im_f = [im_f,f_z];
            subdiff_f = [subdiff_f,subgrad_f(z)];
            subhess_f_cell = [subhess_f_cell(:)',subhess_f(z)];
            
            %%%%%%%%%%%%%%%%%%
            %%%%% Step 3 %%%%%
            %%%%%%%%%%%%%%%%%%
            
            [z,theta,lambda] = solve_subproblem(sample_pts,im_f,subdiff_f,subhess_f_cell,x_arr(:,i),eps,optns);
            f_z = f(z);
            
            %%%%%%%%%%%%%%%%%%
            %%%%% Step 4 %%%%%
            %%%%%%%%%%%%%%%%%%
            if(theta - f_x > -tau*eps)
                reduction_flag = true;
            end
        end
        
        numsample_arr(i) = numsample_arr(i) + N_sample - N_reused;
        
        if(disp_flag > 1)
            disp(['    ',num2str(N_sample - init_N_sample - N_reused),' additional sample points required. ',...
                '(Nonzero multipliers = ',num2str(sum(lambda.ineqnonlin(1:end-1) ~= 0)),...
                '. Bound multiplier = ',num2str(lambda.ineqnonlin(end)),').']);
            disp(['    f(x^i) = ',sprintf('%9.2e', f_x)]);
            disp(['    f(z)   = ',sprintf('%9.2e', f_z)]);
            disp(['    theta  = ',sprintf('%9.2e', theta)]);
            disp(['    f(x^i) - f(z) = ',sprintf('%9.2e', f_x - f_z)]);
            disp(['    min(f(sample_pts)) = ',sprintf('%9.2e', min(im_f))]);
            
            [~,fval] = quadprog(2*(subdiff_f'*subdiff_f),zeros(N_sample,1),[],[],ones(1,N_sample),1,zeros(N_sample,1),ones(N_sample,1),[],optimoptions('quadprog','Display','off'));
            om = sqrt(fval);
            disp(['    Approx. opt. measure = ',sprintf('%9.2e', om)])
        end
            
        if(reduction_flag)
            
            %%%%%%%%%%%%%%%%%%
            %%%%% Step 4 %%%%%
            %%%%%%%%%%%%%%%%%%
            
            if(disp_flag > 1)
                disp(['Reducing eps and tau: ',...
                    '(theta - f(x^i))/eps = ',sprintf('%9.2e', (theta - f_x)/eps),...
                    ]);
                disp(' ');
            end
            eps = kappa_eps*eps;
            tau = kappa_tau*tau;
            eps_ind(j_eps+1) = i;
            j_eps = j_eps+1;
        else
            
            %%%%%%%%%%%%%%%%%%
            %%%%% Step 6 %%%%%
            %%%%%%%%%%%%%%%%%%
            
            if(disp_flag > 1)
                disp(['Point accepted: (theta - f(x))/eps = ',sprintf('%9.2e', (theta - f_x)/eps),...
                    ', ||z - x^i|| = ',sprintf('%9.2e', norm(z - x_arr(:,i),2)),...
                    ])
                disp(' ');
            end

            lambda_arr(:,i) = lambda.ineqnonlin(N_sample+1);
            x_arr(:,i+1) = z;
            f_x = f_z;
            i = i+1;
        end
    end
        
    x_arr = x_arr(:,1:i);
    if(reduction_flag)
        numsample_arr = numsample_arr(1:i);
    else
        numsample_arr = numsample_arr(1:(i-1));
    end
    lambda_arr = lambda_arr(:,1:(i-1));
    eps_ind = eps_ind(1:j_eps);
    fval = f_x;
    
    if(disp_flag > 0)
        disp('Algorithm finished!');
        disp(['    Final objective value            = ',sprintf('%9.2e', fval)]);
        disp(['    Required subgradient evaluations = ',num2str(sum(numsample_arr))]);
    end
    if(i > max_iter)
        disp('Warning: Maximal number of iterations reached. Likely no convergence.')
    end

end

