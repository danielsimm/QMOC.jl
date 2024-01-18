using QMOC
using QuantumClifford
using Statistics
using DelimitedFiles

L = 12
trajectories = [QMOC.trajectory(:YaoKivelsonOrientable, L, 6*L^2, "test", [p, 1-p], false, :debug, 1, 1, 1, 1) for p in [0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0, 0.63]]

inits = [QMOC._DHC_ZZ_operators(L)[1:end-1]..., QMOC._DHC_largeloop_operators(L)[1:end-1]..., QMOC._DHC_smallloop_operators(L)..., QMOC._DHC_wilsonline_operators(L)...]
operators = QMOC.get_operators(trajectories[1])

output = []
Threads.@threads for j in 1:500
    tmis = zeros(8, 3*L+1)
    for i in eachindex(trajectories)
        state = Stabilizer(inits)
        tmis[i, 1] = QMOC.tmi(state, trajectories[i])
        for time in 1:3*L
            QMOC.circuit!(state, trajectories[i], operators)
            tmis[i, time+1] = QMOC.tmi(state, trajectories[i])
        end
    end
    push!(output, tmis)
end

avg = mean(output)
err = std(output) ./ sqrt(length(output))

writedlm("equilibration_times_avg.txt", avg)
writedlm("equilibration_times_err.txt", err)


