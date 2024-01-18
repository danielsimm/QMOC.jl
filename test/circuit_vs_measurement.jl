using QMOC
using BenchmarkTools
using Statistics
using DelimitedFiles

Ls = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48]
trajectories = [QMOC.trajectory(
    :YaoKivelsonNonorientable,
    L,
    6*L^2,
    "benchmark",
    [0.5, 0.5],
    false,
    :none,
    1,
    3*L,
    1,
    1) for L in Ls]

states = QMOC.initialise.(trajectories)
operators = QMOC.get_operators.(trajectories)

### thermalise states ###
Threads.@threads for i in eachindex(trajectories)
    for t in 1:trajectories[i].thermalization_steps
        QMOC.circuit!(states[i], trajectories[i], operators[i])
    end
end

function measurement(state, trajectory)
    meas = Measurement(QMOC.entropy(state, trajectory), QMOC.tmi(state, trajectory))
end


circuit_benchmarks = []
measurement_benchmarks = []

for i in eachindex(trajectories)
    push!(circuit_benchmarks, @benchmark QMOC.circuit!($states[$i], $trajectories[$i], $operators[$i]) seconds=600)
    push!(measurement_benchmarks, @benchmark measurement($states[$i], $trajectories[$i]) seconds=600)
end

circuit_benchmarks
median(measurement_benchmarks[3]).time

# write out benchmark results

median_times = hcat("median_time (ns)", [median(benchmark).time for benchmark in circuit_benchmarks]...)
mean_times = hcat("mean_time (ns)", [mean(benchmark).time for benchmark in circuit_benchmarks]...)
std_times = hcat("std_time (ns)", [std(benchmark).time for benchmark in circuit_benchmarks]...)
min_times = hcat("min_time (ns)", [minimum(benchmark).time for benchmark in circuit_benchmarks]...)
max_times = hcat("max_time (ns)", [maximum(benchmark).time for benchmark in circuit_benchmarks]...)
memory_estimate = hcat("memory_estimate (byte)", [benchmark.memory for benchmark in circuit_benchmarks]...)

result = vcat(median_times, mean_times, std_times, min_times, max_times, memory_estimate)
writedlm("circuit_benchmarks.txt", result)

median_times = hcat("median_time (ns)", [median(benchmark).time for benchmark in measurement_benchmarks]...)
mean_times = hcat("mean_time (ns)", [mean(benchmark).time for benchmark in measurement_benchmarks]...)
std_times = hcat("std_time (ns)", [std(benchmark).time for benchmark in measurement_benchmarks]...)
min_times = hcat("min_time (ns)", [minimum(benchmark).time for benchmark in measurement_benchmarks]...)
max_times = hcat("max_time (ns)", [maximum(benchmark).time for benchmark in measurement_benchmarks]...)
memory_estimate = hcat("memory_estimate (byte)", [benchmark.memory for benchmark in measurement_benchmarks]...)

result = vcat(median_times, mean_times, std_times, min_times, max_times, memory_estimate)
writedlm("measurement_benchmarks.txt", result)