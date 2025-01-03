## Script which applies the VUbundle algorithm from https://github.com/GillesBareilles/NonSmoothSolvers.jl to the half-and-half problem

using NonSmoothSolvers
using NonSmoothProblems
using Plots
using MAT
using Dates

n = 8
pb = Halfhalf()
x = zeros(8) .+ 20.08
optparams = OptimizerParams(iterations_limit=10^5, trace_length=20)

## Run VU bundle method
o = VUbundle()
xfinal_vu, tr = optimize!(pb, o, x; optparams);

## Process output
num_evals = length(tr[end].additionalinfo[8])
xl_arr = zeros(n,num_evals)
fxl_arr = zeros(num_evals)

for i in 1:num_evals
    xl_arr[:,i] = getindex(tr[end].additionalinfo[8],i)
    fxl_arr[i] = sqrt(xl_arr[:,i]'*pb.A*xl_arr[:,i]) + xl_arr[:,i]'*pb.B*xl_arr[:,i]
end

num_iter = length(tr)
fxj_arr = zeros(num_iter)
bbcalls_arr = zeros(num_iter)

for i in 1:num_iter
    fxj_arr[i] = tr[i].Fx;

    if(i > 1)
        bbcalls_arr[i] = tr[i].additionalinfo[4];
    end
end

p1 = plot(1:num_iter,log10.(fxj_arr),markershape=:circle, markersize=2, legend=false);
xlabel!("Iter")
ylabel!("Distance to minimal value")

p2 = plot(bbcalls_arr,log10.(fxj_arr),markershape=:circle, markersize=2, legend=false)
xlabel!("Subgrad. eval.")

p3 = plot(1:num_evals,log10.(fxl_arr),markershape=:circle, markersize=2, legend=false)
xlabel!("Subgrad. eval.")

plot(p1,p2,p3,layout=(1,3))

# Export result to mat-file
file = matopen("julia-results-" * Dates.format(Dates.now(), dateformat"dd-mm-yy_HH-MM-SS") * ".mat", "w")
write(file, "subgrad_arr_vu", bbcalls_arr)
write(file, "f_arr_vu", fxj_arr)
close(file)
