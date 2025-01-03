function [f_x,subgrad_f_x] = my_fun(x)

    global function_data
    f_x = function_data.f(x);
    subgrad_f_x = function_data.subgrad_f(x);  

end

