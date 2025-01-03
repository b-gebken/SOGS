function v = my_slqpgs_problem(x,d,j,o)

% Switch on evaluation option
switch o
  
  case 0 % Objective function value
    
    v = d.f(x);
    
  case 1 % Objective gradient value (as column vector)
    
    v = d.subgrad_f(x);
    
  case 2 % j-th equality constraint value
    
  case 3 % j-th equality constraint gradient value (as row vector)
    
  case 4 % j-th inequality constraint value
    
  case 5 % j-th inequality constraint gradient value (as row vector)

  otherwise
    
    % This error is not supposed to happen!
    error('SLQP-GS: Option not given to problem function evaluator.');
    
end