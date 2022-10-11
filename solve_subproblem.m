function [z,theta,lambda] = solve_subproblem(sample_pts,im_f,subdiff_f,subhess_f_cell,x0,eps,optns)
%   Solves the subproblem (17) from [1] via fmincon. 
%   Compared to the version in the paper, the problem that is solved here 
%   has additional box constraints 
%       z in [x0 - eps, x0 + eps]. 
%   These are redundant from a mathematical point of view, since we already
%   impose 
%        z in B_eps(x0)
%   with respect to the 2-norm, but they led to better behavior of fmincon
%   in our tests.
%
%   [1] Gebken, "Using second-order information in gradient sampling
%   methods for nonsmooth optimization" (2022)

n = size(subdiff_f,1);

fun = @(in) get_fun(in,n);
nonlcon = @(in) get_nonlcon(in,n,sample_pts,im_f,subdiff_f,subhess_f_cell,x0,eps);
laghess = @(in,lambda) get_laghess(in,lambda,n,subhess_f_cell);

optns.SpecifyObjectiveGradient = true;
optns.SpecifyConstraintGradient = true;

if(strcmp(optns.Algorithm,'interior-point'))
    optns.HessianFcn = laghess;
end

exitflag = -Inf;
rand_flag = false;
attempt_counter = 0;
while(exitflag <= 0)
    if(rand_flag)
        rnd_dir = 2*rand(n,1)-1; rnd_dir = rnd_dir/norm(rnd_dir,2);
        x_rand = x0 + 0.5*eps*rnd_dir;
        [sol,~,exitflag,output,lambda] = fmincon(fun,[x_rand;max(nonlcon([x_rand;0])) + 1],[],[],[],[],[(x0 - eps);-Inf],[(x0 + eps);Inf],nonlcon,optns);
    else
        [sol,~,exitflag,output,lambda] = fmincon(fun,[x0;max(nonlcon([x0;0])) + 1],[],[],[],[],[(x0 - eps);-Inf],[(x0 + eps);Inf],nonlcon,optns);
    end
    attempt_counter = attempt_counter + 1;
    
    if(exitflag == 0 && output.constrviolation < optns.ConstraintTolerance)
        disp('Warning: Feasible point with exitflag = 0.')
        break
    end
    
    if(exitflag <= 0)
        disp(['Warning: Subproblem solution failed (attempt ',num2str(attempt_counter),'), retrying with random starting point (flag = ',num2str(exitflag),', ||z|| = ',num2str(norm(sol(1:n) - x0,2)),', eps = ',num2str(eps),').'])
        rand_flag = true;
    end
end

z = sol(1:n);

% Since the inequality constraints corresponding to the second-order Taylor
% expansions are not exactly 0 when they are active (due to numerical
% inaccuracies), we can get a better objective value than 
%   fun(sol) = sol(n+1)
% by taking the actual maximum of the Taylor expansions "by hand":
k = numel(subhess_f_cell);
nonlcon_sol = nonlcon(sol);
theta = max(nonlcon_sol(1:k) + sol(n+1));

end

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