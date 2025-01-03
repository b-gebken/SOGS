% A wrapper for recording all point in which the subgradient was evaluated

function in = eval_pts_wrapper(in)

    global eval_pts_global
    eval_pts_global(:,end+1) = in;

end

