struct Simulation
    type::Symbol
    size::Int
    nqubits::Int
    name::String
    checkpoints::Bool
    verbosity::Symbol
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
    num_trajectrories::Int64
    parameter_set::Vector{Vector{Real}}
    ensemble::Vector{Vector{Trajectory}}
end

function simulation(
    type::Symbol,
    size::Int64,
    name::String,
    num_trajectories::Int64,
    parameter_set;
    checkpoints::Bool=false,
    verbosity::Symbol=:info,
    nqubits::Int64=_number_of_qubits(type, size),
    thermalization_steps::Int64=3 * size,
    number_of_measurements::Int64=100,
    measurement_steps::Int64=3)
    ensemble = Vector{Vector{Trajectory}}(undef, length(parameter_set))
    for i in eachindex(ensemble)
        ensemble[i] = Vector{Trajectory}(undef, num_trajectories)
        for j in eachindex(ensemble[i])
            ensemble[i][j] = trajectory(type, size, nqubits, name, parameter_set[i], checkpoints, verbosity, j, thermalization_steps, measurement_steps, number_of_measurements)
        end
    end
    return Simulation(type, size, nqubits, name, checkpoints, verbosity, thermalization_steps, measurement_steps, number_of_measurements, num_trajectories, parameter_set, ensemble)
end

function test_simulation(type)
    return simulation(type, 12, "test_$(type)", 4, [[1//3, 1//3, 1//3], [1//2, 1//2, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]; checkpoints=true, verbosity=:debug)
end

function test_simulation_large(type)
    return simulation(type, 12, "test_$(type)_large", 10, parameter_line(:center, :px, 15); checkpoints=true)
end

function trajectory(
    type,
    size,
    nqubits,
    name,
    parameters,
    checkpoints,
    verbosity,
    index,
    thermalization_steps,
    measurement_steps,
    number_of_measurements) ::Trajectory
    if type == :ChainPP
        return PPChainTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :ChainPQ
        return PQChainTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :ChainKekule
        return KekuleChainTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :Kekule
        return KekuleTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :Kitaev
        return KitaevTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :YaoKivelsonXYZ
        return YaoKivelsonXYZTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :YaoKivelsonNonorientable
        return YaoKivelsonNonorientableTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :YaoKivelsonOrientable
        return YaoKivelsonOrientableTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    else
        error("Type $(type) not implemented.")
    end
end

function test_trajectory(type)
    return trajectory(type, 12, _number_of_qubits(type, 12), "test_$(type)", [1//3, 1//3, 1//3], true, :debug, 1, 36, 3, 10)
end

function _number_of_qubits(type::Symbol, size::Int) ::Int
    if type in [:ChainPP, :ChainPQ, :ChainKekule]
        return size
    elseif type in [:Kekule, :Kitaev]
        return 2 * size^2
    elseif type in [:YaoKivelsonXYZ, :YaoKivelsonNonorientable, :YaoKivelsonOrientable]
        return 3 * size^2
    else
        error("Type $(type) not implemented.")
    end
end

# function simulate(simulation::Simulation)

#     # set BLAS threads to 1 to avoid oversubscription
#     BLAS.set_num_threads(1)
#     trajectories = []
#     for i in eachindex(simulation.ensemble)
#         for j in eachindex(simulation.ensemble[i])
#             push!(trajectories, simulation.ensemble[i][j])
#         end
#     end

#     pmap(run, trajectories; retry_delays = zeros(3))

#     # commit metadata to file
#     commitMetadata(simulation)

# end

function simulate(sim::Simulation)
    
    MPI.Init()


    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    world_size = MPI.Comm_size(comm)
    nworkers = world_size - 1

    root = 0

    # collect all trajectories
    trajectories = []
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            push!(trajectories, sim.ensemble[i][j])
        end
    end
    trajectories = unique(trajectories)
    ntrajectories = length(trajectories)

    # skip archived trajectories
    if rank == root
        archived_hashes = archivedTrajectories()
        for i in eachindex(trajectories)
            if hash(trajectories[i]) in archived_hashes
               trajectories[i] = nothing
            end
        end

    end

    MPI.Barrier(comm)
    if rank == root
        println("Starting $(sim.name) with $(ntrajectories) trajectories on $(nworkers) workers...")
    end
    
    # distribute indices to workers
    if rank == root
        # randomly distribute indices
        indices = randperm(ntrajectories)
        part = [indices[i:nworkers:end] for i in 1:nworkers]
        for i in 1:nworkers
            MPI.send(part[i],comm; dest=i)
        end
    else
        todo = MPI.recv(comm)
        for i in todo
            run(trajectories[i])
        end
    end
    if rank == root
        println("Dispatched to $(nworkers) MPI procs. Waiting for results...")
    end
    MPI.Barrier(comm)
    if rank == root
        if boolComplete(sim)
            println("All trajectories completed.")
            writeMetadata(sim)
            println("Finished $(sim.name) with $(ntrajectories) trajectories on $(nworkers) workers.")
        else
            println("Not all trajectories completed.")
        end
    end
    
    MPI.Barrier(comm)
    MPI.Finalize()
end

Base.show(io::IO, sim::Simulation) = 
print(io,
"--- ", sim.name, " ---", "\n",
"Type: ", sim.type, "\n",
"Size: ", sim.size, "\n",
"Number of qubits: ", sim.nqubits, "\n",
"Checkpoints: ", sim.checkpoints, "\n",
"Verbosity: ", sim.verbosity, "\n",
"Thermalization steps: ", sim.thermalization_steps, "\n",
"Measurements: ", sim.number_of_measurements, " @ ", sim.measurement_steps, " steps", "\n",
"Number of trajectories: ", sim.num_trajectrories, "\n")

Base.show(io::IO, traj::Trajectory) =
print(io,
"--- ", typeof(traj), " ---", "\n",
"Name: ", traj.name, "\n",
"Size: ", traj.size, "\n",
"Number of qubits: ", traj.nqubits, "\n",
"Checkpoints: ", traj.checkpoints, "\n",
"Verbosity: ", traj.verbosity, "\n",
"Index: ", traj.index, "\n",
"Thermalization steps: ", traj.thermalization_steps, "\n",
"Measurements: ", traj.number_of_measurements, " @ ", traj.measurement_steps, " steps", "\n",
"Parameters: ", traj.params, "\n")