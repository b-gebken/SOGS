% Adds jet elements to the memory-struct. If the maximum size of the memory
% is reached, then the oldest elements are removed. 

function memory = add_to_memory(new_sample_pts,new_f_vals,new_subgrads,new_subhess,memory)

    memory.sample_pts = [memory.sample_pts,new_sample_pts];
    memory.f_vals = [memory.f_vals,new_f_vals];
    memory.subgrads = [memory.subgrads,new_subgrads];
    memory.subhess = [memory.subhess,new_subhess];
    
    if(size(memory.sample_pts,2) > memory.max_size)
        memory.sample_pts = memory.sample_pts(:,end-memory.max_size+1:end);
        memory.f_vals = memory.f_vals(end-memory.max_size+1:end);
        memory.subgrads = memory.subgrads(:,end-memory.max_size+1:end);
        memory.subhess = memory.subhess(end-memory.max_size+1:end);
    end
end

