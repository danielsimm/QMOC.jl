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
    verbosity::Symbol=:low,
    nqubits::Int64 = _number_of_qubits(type, size),
    thermalization_steps::Int64=3*size,
    number_of_measurements::Int64=100,
    measurement_steps::Int64=3
    )
    ensemble = Vector{Vector{Trajectory}}(undef, length(parameter_set))
    for i in eachindex(ensemble)
        ensemble[i] = Vector{Trajectory}(undef, num_trajectories)
        for j in eachindex(ensemble[i])
            ensemble[i][j] = trajectory(type, size, nqubits, name, parameter_set[i], checkpoints, verbosity, j, thermalization_steps, measurement_steps, number_of_measurements)
        end
    end
    return Simulation(type, size, nqubits, name, checkpoints, verbosity, thermalization_steps, measurement_steps, number_of_measurements, num_trajectories, parameter_set, ensemble)
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
    number_of_measurements)
    if type == :ChainPP
        return PPChainTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :ChainPQ
        return PQChainTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :ChainPPFast
        return PPChainTrajectoryFast(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :Kekule
        return KekuleTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :Kitaev
        return KitaevTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :YaoKivelsonXYZ
        return YaoKivelsonXYZTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    elseif type == :YaoKivelsonJJ
        return YaoKivelsonJJTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    else
        error("Type $(type) not implemented.")
    end
end

function _number_of_qubits(type::Symbol, size::Int)
    if type in [:ChainPP, :ChainPQ, :ChainPPFast]
        return size
    elseif type in [:Kekule, :Kitaev]
        return 2*size^2
    elseif type in [:YaoKivelsonXYZ, :YaoKivelsonJJ]
        return 3*size^2
    else
        error("Type $(type) not implemented.")
    end
end


function simulate(simulation::Simulation)
    
    # set BLAS threads to 1 to avoid oversubscription
    BLAS.set_num_threads(1)
    trajectories = []
    for i in eachindex(simulation.ensemble)
        for j in eachindex(simulation.ensemble[i])
            push!(trajectories, simulation.ensemble[i][j])
        end
    end
    
    pmap(run, trajectories)
end

