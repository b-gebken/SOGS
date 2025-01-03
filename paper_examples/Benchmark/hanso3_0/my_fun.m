function [f_x,subgrad_f_x] = my_fun(x,pars)

    f_x = pars.my_f(x);
    if(nargout > 1)
        subgrad_f_x = pars.my_subgrad_f(x);  
    end

end

