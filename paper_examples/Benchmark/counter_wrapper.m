% A wrapper for counting the number of evaluations

function in = counter_wrapper(in,counter_flag)

    if(counter_flag == 1)
        global f_counter
        f_counter = f_counter + 1;
    elseif(counter_flag == 2)
        global subgrad_counter
        subgrad_counter = subgrad_counter + 1;
    else
        global subhess_counter
        subhess_counter = subhess_counter + 1;
    end

end

