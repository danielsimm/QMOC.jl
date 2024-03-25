using QMOC
using DelimitedFiles
Ls = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48]
trajectories = [QMOC.trajectory(:YaoKivelsonNonorientable, L, 6*L^2, "test", [0.5, 0.5], false, :none, 1, L, 0, 0) for L in [12, 16, 20, 24, 28, 32, 36, 40, 44, 48]]

times = zeros(length(trajectories))
for _ in 1:5
    for (i, trajectory) in enumerate(trajectories)
        isfile("data/measurements/$(hash(trajectories[i])).jld2") && rm("data/measurements/$(hash(trajectories[i])).jld2")
        start = time()
        run(trajectory)
        stop = time()
        times[i] += (stop - start)
        isfile("data/measurements/$(hash(trajectories[i])).jld2") && rm("data/measurements/$(hash(trajectories[i])).jld2")
    end
end

times = times ./5
to_file = hcat(Ls, times)
writedlm("data/performance_bench2.txt", to_file)


