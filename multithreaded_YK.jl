using QMOC
using DelimitedFiles
using QuantumClifford
using MPI
using Dates
using Random

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

parameter_set = parameter_set_orientable
mode = :YaoKivelsonOrientable
n_samples = Threads.nthreads()

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

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
world_size = MPI.Comm_size(comm)
if rank == 0
    indices = randperm(length(parameter_set))
    part = [indices[i:(world_size):end] for i in 1:(world_size)]
    for i in 1:world_size-1
        MPI.send(part[i+1],comm; dest=i)
    end
    todo = part[1]
    println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- starting...")
    for i in todo
        sample_trajectory(L, mode, parameter_set[i], n_samples, L+10, 100, "scratchdata/YK_$(mode)_$(L)_$(i)")
    end
    println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- done.")
else
    todo = MPI.recv(comm)
    println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- starting...")
    for i in todo
        sample_trajectory(L, mode, parameter_set[i], n_samples, L+10, 100, "scratchdata/YK_$(mode)_$(L)_$(i)")
    end
    println("$(Dates.format(now(), "HH:MM:SS")) -- $(mode) | $(L) | idx $(todo) | rank $(rank) -- done.")
end
MPI.Barrier(comm)
MPI.Finalize()