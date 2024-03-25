using LinearAlgebra
BLAS.set_num_threads(1)

using QMOC
using DelimitedFiles
using QuantumClifford

using Dates
using Random


L = 4
function K_range(L, nu, Kc, ΔK)
    K_min = max(0, Kc - ΔK/(L^(1/nu)))
    K_max = Kc + ΔK/(L^(1/nu))
    return LinRange(K_min, K_max, 50)
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

parameter_set = parameter_set_nonorientable
mode = :YaoKivelsonNonorientable
n_samples = Threads.nthreads()

trajectory = QMOC.trajectory(mode, L, 6*L^2, "test", parameter_set[1], false, :debug, 1, 1, 1, 1)
state = QMOC.initialise(trajectory)
operators = QMOC.get_operators(trajectory)
this_state = copy(state)
for t in 1:2*L
    QMOC.circuit!(this_state, trajectory, operators)
end
tab = this_state.tab
for i in eachindex(tab)
    println(tab[i])
end
for t in 1:subsamples
    QMOC.circuit!(this_state, trajectory, operators)
    tmi = QMOC.tmi(this_state, trajectory)
    open("$(outputname)_sub$(s).tmi", "a+") do io
        writedlm(io, tmi)
    end
end

function sample_trajectory(L, mode, params, samples, thermalization_time, subsamples, outputname)
    trajectory = QMOC.trajectory(mode, L, 6*L^2, outputname, params, false, :debug, 1, 1, 1, 1)
    state = QMOC.initialise(trajectory)
    operators = QMOC.get_operators(trajectory)

    Threads.@threads for s in 1:samples
        this_state = copy(state)
        for t in 1:thermalization_time
            QMOC.circuit!(this_state, trajectory, operators)
        end
        for t in 1:subsamples
            QMOC.circuit!(this_state, trajectory, operators)
            tmi = QMOC.tmi(this_state, trajectory)
            open("$(outputname)_sub$(s).tmi", "a+") do io
                writedlm(io, tmi)
            end
        end
    end 
end



println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- starting...")
sample_trajectory(L, mode, parameter_set[ind], n_samples, L+15, 200, "data_nonorientable_64/$(ind)/$(rank+1)")
println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- done.")

MPI.Barrier(comm)
MPI.Finalize()