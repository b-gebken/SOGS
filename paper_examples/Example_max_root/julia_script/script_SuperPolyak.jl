## Script which applies SuperPolyak from https://github.com/COR-OPT/SuperPolyak.jl to the function (5.2)

using LinearAlgebra
using SuperPolyak
using Plots
using MAT
using Dates

include("max_root.jl")

n = 100
f = max_root
subgrad = grad_max_root
x0 = max_root_x0(n)

result = SuperPolyak.superpolyak(
    f,
    subgrad,
    x0,
)

## Results
x_min = zeros(n)
f_min = f(x_min)

num_iter = size(result.loss_history,1)-1
println("||x - x^*|| = ",norm(result.solution - x_min,2))
println("f(x)        = ",f(result.solution))
println("#iter       = ",num_iter)
println("#subgrad    = ",sum(result.oracle_calls))

## Plots
p1 = plot(1:num_iter+1,log10.(result.loss_history .- f_min),markershape=:circle, markersize=2)
p2 = plot(1:num_iter+1,result.oracle_calls,markershape=:circle, markersize=2)
p3 = plot(tril(ones(num_iter+1,num_iter+1))*result.oracle_calls, log10.(result.loss_history .- f_min), markershape=:circle, markersize=2)
plot(p1,p2,p3,layout=(1,3))

# Export result to mat-file
file = matopen("julia-results-" * Dates.format(Dates.now(), dateformat"dd-mm-yy_HH-MM-SS") * ".mat", "w")
write(file, "subgrad_arr_SupPol", result.oracle_calls)
write(file, "f_arr_SupPol", result.loss_history)
close(file)