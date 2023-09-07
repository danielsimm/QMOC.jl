using LinearAlgebra
using Statistics
using BenchmarkTools
include("LatticeCircuits.jl")

# mkpath("data")
# touch("data/benchmarks.jld2")


### Measurement ###

function HC_timestep_suite(length::Int64, step=1)
    Ls = [i for i in 1:step:length]
    # check for existing data

    try
        jldopen("data/benchmarks.jld2", "r") do file
            if haskey(file, "HC_timestep_times")
                times = file["HC_timestep_times"]
                if length(times) == length(Ls)
                    return Ls, times
                end
            end
        end
    catch

    end
    


    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1])
    scatterlines!(ax, Ls, Float64.(times), color=:blue)
    fit_res = power_fit(Ls[3:end], Float64.(times)[3:end])
    a = fit_res[1]
    b = fit_res[2]
    lines!(ax, 1:40, a.*(1:40).^b, color=:red)
    fig
end


function benchmark_HC_timestep(L::Int64) # returns the time (seconds) and memory usage (KiB) of a single timestep
    init = HC_initial_state(L)
    for i in 1:L
        kitaev_timestep!(init, 1//3, 1//3, 1//3)
    end
    result = @benchmark kitaev_timestep!($init, 1//3, 1//3, 1//3)
    return mean(result).time / 1e9, mean(result).memory / 1e6
end

# Ls = [i for i in 1:20]
# times = []
# for i in eachindex(Ls)
#     push!(times, benchmark_HC_timestep(Ls[i])[1])
# end

# times = times .* 1e3 # convert to milliseconds

# using CairoMakie
# using CurveFit

function BLAS_benchmark(L)
    init = HC_initial_state(L)
    for i in 1:3*L
        kitaev_timestep!(init, 1//3, 1//3, 1//3)
    end
    times_evolution = []
    times_TMI = []
    state = deepcopy(init)
    BLAS.set_num_threads(1)
    push!(times_evolution, @benchmark kitaev_timestep!($state, 1//3, 1//3, 1//3))
    push!(times_TMI, @benchmark HC_TMI($state, $L))

    state = deepcopy(init)
    BLAS.set_num_threads(2)
    push!(times_evolution, @benchmark kitaev_timestep!($state, 1//3, 1//3, 1//3))
    push!(times_TMI, @benchmark HC_TMI($state, $L))

    state = deepcopy(init)
    BLAS.set_num_threads(3)
    push!(times_evolution, @benchmark kitaev_timestep!($state, 1//3, 1//3, 1//3))
    push!(times_TMI, @benchmark HC_TMI($state, $L))

    state = deepcopy(init)
    BLAS.set_num_threads(4)
    push!(times_evolution, @benchmark kitaev_timestep!($state, 1//3, 1//3, 1//3))
    push!(times_TMI, @benchmark HC_TMI($state, $L))

    return times_evolution, times_TMI
end
    
function fast_projection_benchmark(L)
    traj1 = PPChainTrajectory(L, L, "test", [1//3, 1//3, 1//3], false, :off, 1, 100, 1, 100)
    traj2 = PPChainTrajectoryFast(L, L, "test", [1//3, 1//3, 1//3], false, :off, 1, 100, 1, 100)
    init = initialise(traj1)
    for i in 1:3*L
        circuit!(init, traj1)
    end
    state1 = deepcopy(init)
    state2 = deepcopy(init)
    b1 = @benchmark circuit!($state1, $traj1)
    b2 = @benchmark circuit!($state2, $traj2)
    return b1, b2
end
b1, b2 = fast_projection_benchmark(1024)