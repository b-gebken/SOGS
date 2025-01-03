% Function for solving the subproblem (4.6) via IPOPT. For details, see 
%   https://github.com/ebertolazzi/mexIPOPT
%   https://ethz.ch/content/dam/ethz/special-interest/mavt/dynamic-systems-n-control/idsc-dam/Research_Onder/Downloads/IPOPT/IPOPT_MatlabInterface_V0p1.pdf

function [z,theta,mu] = solve_subproblem_IPOPT(W_sample_pts,W_f_vals,W_subgrads,W_subhess,x0,eps,sp_solver_options)

n = size(W_subgrads,1);
k = numel(W_subhess);

% Define funcs struct for IPOPT
funcs.objective = @objective;
funcs.constraints = @constraints;
funcs.gradient = @gradient;
funcs.jacobian = @jacobian;
funcs.jacobianstructure = @jacobianstructure;
funcs.hessian = @hessian;
funcs.hessianstructure = @hessianstructure;

% Define auxdata struct containing the data for the subproblem
auxdata.n = n;
auxdata.k = k;
auxdata.eps = eps;
auxdata.x0 = x0;
auxdata.im_f = W_f_vals;
auxdata.subdiff_f = W_subgrads;
auxdata.sample_x = W_sample_pts; 
auxdata.subhess_f_cell = W_subhess;

% Set IPOPT options
IPOPT_opts.lb = [(x0 - eps);-Inf]';
IPOPT_opts.ub = [(x0 + eps);Inf]';
IPOPT_opts.cl = -Inf(1,k+1);
IPOPT_opts.cu = zeros(1,k+1);

IPOPT_opts.auxdata = auxdata;
IPOPT_opts.ipopt.max_iter = 1000;
IPOPT_opts.ipopt.tol = sp_solver_options.tol;
IPOPT_opts.ipopt.print_level = 0;  

IPOPT_opts.ipopt.line_search_method = 'cg-penalty'; % Solved status problems in 7 and 20

% Define initial point for IPOPT. (x0 is modified to avoid an error in the
% restoration phase of IPOPT, which was likely caused by the gradient of
% the eps-ball constraint being zero in x0.)  
x0_mod = x0 + 1/2*eps*[1;zeros(n-1,1)];
init_pt = [x0_mod',max(funcs.constraints([x0_mod',0],auxdata)) + 1];

% Run the solver
[sol_IPOPT, info_IPOPT] = ipopt_auxdata(init_pt, funcs, IPOPT_opts);

% Check IPOPT status
if(info_IPOPT.status ~= 0)
    fprintf(2,['Warning: IPOPT status error (',num2str(info_IPOPT.status),').\n'])
end

% Prepare the output
z = sol_IPOPT(1:n)';
nonlcon_sol = funcs.constraints(sol_IPOPT,auxdata);
theta = max(nonlcon_sol(1:k) + sol_IPOPT(n+1));
mu.ineqnonlin = info_IPOPT.lambda;

end

% Auxiliary functions
function f = objective(in,auxdata)
    f = in(auxdata.n + 1);
end

function g = gradient(in,auxdata)
    g = [zeros(auxdata.n,1);1]';
end

function c = constraints(in,auxdata)
	k = auxdata.k;
    n = auxdata.n;
    
    tmp1 = in(1:n)' - auxdata.sample_x;
    c = [(auxdata.im_f + sum(auxdata.subdiff_f .* tmp1,1) - in(n+1))'; ...
        0];

    tmp2 = zeros(k,1);
    for i = 1:k
        tmp2(i) = tmp1(:,i)'*auxdata.subhess_f_cell{i}*tmp1(:,i);
    end
    c(1:k) = c(1:k) + 1/2*tmp2;
    c(k+1) = norm(in(1:n)' - auxdata.x0,2)^2 - auxdata.eps^2;

    %%% Slower, but more readable version
    % k = auxdata.k;
    % c = zeros(k+1,1);
    % n = auxdata.n;
    % 
    % for i = 1:k
    %     c(i) = auxdata.im_f(i) + auxdata.subdiff_f(:,i)'*(in(1:n)' - auxdata.sample_x(:,i)) + 1/2*(in(1:n)' - auxdata.sample_x(:,i))'*auxdata.subhess_f_cell{i}*(in(1:n)' - auxdata.sample_x(:,i)) - in(n+1);
    % end
    % c(k+1) = norm(in(1:n)' - auxdata.x0,2)^2 - auxdata.eps^2;
end

function J = jacobian(in,auxdata)
    k = auxdata.k;
    n = auxdata.n;
    
    J = [auxdata.subdiff_f',-ones(k,1);
        2*(in(1:n) - auxdata.x0'),0];
    tmp1 = in(1:n)' - auxdata.sample_x;

    tmp2 = zeros(k,n);
    for i = 1:k
        tmp2(i,:) = tmp1(:,i)'*auxdata.subhess_f_cell{i};
    end
    J(1:k,1:n) = J(1:k,1:n) + tmp2; 
    
    J = sparse(J);

    %%% Slower, but more readable version
    % k = auxdata.k;
    % n = auxdata.n;
    % J = zeros(k+1,n+1);
    % for i = 1:k
    %     J(i,:) = [auxdata.subdiff_f(:,i) + auxdata.subhess_f_cell{i}*(in(1:n)' - auxdata.sample_x(:,i)); -1]';
    % end
    % J(k+1,:) = [2*(in(1:n)' - auxdata.x0);0]';
    % 
    % J = sparse(J);
end

function J_str = jacobianstructure(auxdata)
    k = numel(auxdata.subhess_f_cell);
    n = auxdata.n;
    J_str = sparse(ones(k+1,n+1));
end

function H = hessian (in, sigma, lambda, auxdata)
    k = numel(auxdata.subhess_f_cell);
    n = auxdata.n;
    H = zeros(n+1);

    tmp = zeros(n);
    for i = 1:k
        tmp = tmp + lambda(i)*auxdata.subhess_f_cell{i};
    end
    H = [tmp,zeros(n,1);zeros(1,n+1)];
    
    H(1:n,1:n) = H(1:n,1:n) + lambda(k+1)*2*eye(n);

    H = sparse(tril(H));
end

function H_str = hessianstructure(auxdata)
    n = auxdata.n;
    H_str = sparse([ones(n,n),zeros(n,1); zeros(1,n), 1]);
    H_str = sparse(tril(H_str));
end