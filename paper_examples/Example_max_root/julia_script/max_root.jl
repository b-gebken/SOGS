function max_root(x::Vector{Float64})
    a = 0.1
    return maximum(sqrt.(abs.(x) .+ a) .- sqrt(a))
end

function grad_max_root(x::Vector{Float64})
    a = 0.1
    (~,I) = findmax(sqrt.(abs.(x) .+ a) .- sqrt(a)) 

    n = size(x,1)
    grad = zeros(n);
    
    s = sign(x[I])
    if(s == 0)
        s = 1
    end
    grad[I] = 1/2*(abs(x[I]) + a)^(-1/2) * s

    return grad
end

function max_root_x0(n::Int64)
    return 5*ones(n)
end
