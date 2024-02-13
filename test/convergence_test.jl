using QMOC
using QuantumClifford
using Statistics
using DelimitedFiles

L = 12
function K_range(L, nu, Kc, ΔK)
    K_min = max(0, Kc - ΔK/(L^(1/nu)))
    K_max = Kc + ΔK/(L^(1/nu))
    return LinRange(K_min, K_max, 10)
end
Kc_orientable = 0.631
nu_orientable = 1.0
Kc_nonorientable = 0.654
nu_nonorientable = 0.94
ΔK = 4
Ks_orientable = K_range(L, nu_orientable, Kc_orientable, ΔK)
Ks_nonorientable = K_range(L, nu_nonorientable, Kc_nonorientable, ΔK)
parameter_set_orientable = [[1-K, K] for K in Ks_orientable]
parameter_set_nonorientable = [[1-K, K] for K in Ks_nonorientable]

function TMI_evolution(L, parameter_set, samples, mode)
    trajectories = [QMOC.trajectory(mode, L, 6*L^2, "test", params, false, :debug, 1, 1, 1, 1) for params in parameter_set]
    tmis = zeros(length(trajectories), 3*L+1)
    errors = zeros(length(trajectories), 3*L+1)
    state = QMOC.initialise(trajectories[1])
    operators = QMOC.get_operators(trajectories[1])
    for i in eachindex(trajectories)
        traj = trajectories[i]
        this_tmis = zeros(3*L+1, samples)
        Threads.@threads for j in 1:samples
            this_state = copy(state)
            this_tmis[1, j] = QMOC.tmi(this_state, traj)
            for t in 1:3*L
                QMOC.circuit!(this_state, traj, operators)
                this_tmis[t+1, j] = QMOC.tmi(this_state, traj)
            end
        end
        tmis[i, :] = mean(this_tmis, dims=2)
        errors[i, :] = std(this_tmis, dims=2) ./ sqrt(samples)
    end
    writedlm("tmis_$(mode)_$(L).txt", tmis)
    writedlm("errors_$(mode)_$(L).txt", errors)
end


TMI_evolution(L, parameter_set_orientable, 100, :YaoKivelsonOrientable)

