% Function for solving the subproblem (4.6) via fmincon. (Less optimized
% than the method using the IPOPT solver.) 

function [z,theta,mu] = solve_subproblem_fmincon(W_sample_pts,W_f_vals,W_subgrads,W_subhess,x0,eps,optns) 

n = size(W_subgrads,1);

% Define function handles for fmincon
fun = @(in) get_fun(in,n);
nonlcon = @(in) get_nonlcon(in,n,W_sample_pts,W_f_vals,W_subgrads,W_subhess,x0,eps);
laghess = @(in,lambda) get_laghess(in,lambda,n,W_subhess);

% Set options for fmincon
optns.SpecifyObjectiveGradient = true;
optns.SpecifyConstraintGradient = true;
if(strcmp(optns.Algorithm,'interior-point'))
    optns.HessianFcn = laghess;
end

% Run the solver
[sol,~,exitflag,~,mu] = fmincon(fun,[x0;max(nonlcon([x0;0])) + 1],[],[],[],[],[(x0 - eps);-Inf],[(x0 + eps);Inf],nonlcon,optns);

% Check the exitflag
if(exitflag < 1)
    fprintf(2,['Warning: fmincon exitflag = ',num2str(exitflag),'.\n'])
end

% Prepare output
z = sol(1:n);
k = numel(W_subhess);
nonlcon_sol = nonlcon(sol);
theta = max(nonlcon_sol(1:k) + sol(n+1));

end

% Auxiliary functions
function [y,grady] = get_fun(in,n)
    y = in(n+1);
    grady = [zeros(n,1);1];
end

function [c,ceq,gradc,gradceq] = get_nonlcon(in,n,sample_x,im_f,subdiff_f,subhess_f_cell,x0,eps)
    ceq = 0;
    gradceq = zeros(n+1,1);
    
    k = numel(subhess_f_cell);
    c = zeros(k+1,1);
    gradc = zeros(n+1,k+1);
    for i = 1:k
        c(i) = im_f(i) + subdiff_f(:,i)'*(in(1:n) - sample_x(:,i)) + 1/2*(in(1:n) - sample_x(:,i))'*subhess_f_cell{i}*(in(1:n) - sample_x(:,i)) - in(n+1);

        gradc(:,i) = [subdiff_f(:,i) + subhess_f_cell{i}*(in(1:n) - sample_x(:,i)); -1];
    end
    c(k+1) = norm(in(1:n) - x0,2)^2 - eps^2;
    gradc(:,k+1) = [2*(in(1:n) - x0);0];
end

function laghess = get_laghess(in,lambda,n,subhess_f_cell)
    
    laghess = zeros(n+1);
    k = numel(subhess_f_cell);
    for i = 1:k
        laghess = laghess + lambda.ineqnonlin(i)*[subhess_f_cell{i},zeros(n,1);zeros(1,n+1)];
    end
    laghess = laghess + lambda.ineqnonlin(k+1)*[2*eye(n),zeros(n,1); zeros(1,n),0];
end