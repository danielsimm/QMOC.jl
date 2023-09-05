abstract type Trajectory end
abstract type HoneycombTrajectory <: Trajectory end
abstract type DecoratedHoneycombTrajectory <: Trajectory end
struct KitaevTrajectory <: HoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

struct KekuleTrajectory <: HoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

struct YaoKivelsonXYZTrajectory <: DecoratedHoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

struct YaoKivelsonJJTrajectory <: DecoratedHoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end




"""
    trajectory(id, init, L, parameters, mode, checkpoints, verbose)

TBW
"""

# function _get_checkpoint(id)
#     try
#         load_dict = load("data/checkpoints/$(id).jld2")
#         return load_dict["state"]
#     catch
#         @warn "Checkpoint $(id) could not be loaded."
#         return nothing
#     end
# end

# function _set_checkpoint(id, state)
#     jldopen("data/checkpoints/$(id).jld2", "w") do file
#         file["state"] = state
#     end
# end

function run(trajectory::Trajectory)
    time = 0
    # initialise
    state = initialise(trajectory)
    # thermalize
    thermalize!(state, trajectory, time)
    # measure
    measure!(state, trajectory, time)
end







    


function _yaokivelson_XYZ_trajectory(id, init, L, parameters, mode, checkpoints, verbose)
    px = parameters[1]
    py = parameters[2]
    pz = parameters[3]
    measurements = 100


    ### Measurements
    entropy = zeros(measurements, L+1)
    TMI = zeros(measurements)
    for measurement in 1:measurements
        for layer in 1:3
            _DHC_XYZ_timestep!(state, px, py, pz)
        end
        Threads.@spawn entropy[measurement, :] = DHC_EE(copy(state), L)
        Threads.@spawn TMI[measurement] = DHC_TMI(copy(state), L)
    end

    return sum(entropy, dims=1)./measurements, sum(TMI)./measurements
end

function _yaokivelson_JJ_trajectory(id, init, L, parameters, mode, checkpoints, verbose)
    pJ = parameters[1]
    measurements = 100

    # Check for existing checkpoints, thermalisation
    if checkpoints && isfile("data/checkpoints/$(id).jld2")
        state = _get_checkpoint(id)
        if state === nothing
            @goto loadErrThermalize
        else
            verbose ? (@info "Successfully loaded checkpoint $(id).") : nothing
        end
    else
        if checkpoints
            verbose ? (@info "No checkpoint found for $(id). Starting thermalisation.") : nothing
        end
        @label loadErrThermalize
        state = init
        for layer in 1:6*L
            _DHC_JJ_timestep!(state, pJ)
        end
        verbose ? (@info "Thermalisation complete.") : nothing
        if checkpoints
            _set_checkpoint(id, state)
            verbose ? (@info "Checkpoint $(id) saved.") : nothing
        end
    end

    ### Measurements
    entropy = zeros(measurements, L+1)
    TMI = zeros(measurements)
    for measurement in 1:measurements
        for layer in 1:3
            _DHC_JJ_timestep!(state, pJ)
        end
        entropy[measurement, :] = DHC_EE(state, L)
        TMI[measurement] = DHC_TMI(state, L)
    end

    return sum(entropy, dims=1)./measurements, sum(TMI)./measurements
end

function _kitaev_trajectory(id, init, L, parameters, mode, checkpoints, verbose)
    
    # Setup
    if mode == :Kitaev
        circuit! = kitaev_timestep!
    elseif mode == :Kekule
        circuit! = kekule_timestep!
    end
    px = parameters[1]
    py = parameters[2]
    pz = parameters[3]
    measurements = 100


    # Check for existing checkpoints, thermalisation
    if checkpoints && isfile("data/checkpoints/$(id).jld2")
        state = _get_checkpoint(id)
        if state === nothing
            @goto loadErrThermalize
        else
            verbose ? (@info "Successfully loaded checkpoint $(id).") : nothing
        end
    else
        if checkpoints
            verbose ? (@info "No checkpoint found for $(id). Starting thermalisation.") : nothing
        end
        @label loadErrThermalize
        state = init
        for layer in 1:3*L
            circuit!(state, px, py, pz)
        end
        verbose ? (@info "Thermalisation complete.") : nothing
        if checkpoints
            _set_checkpoint(id, state)
            verbose ? (@info "Checkpoint $(id) saved.") : nothing
        end
    end

    ### Measurements
    entropy = zeros(measurements, L+1)
    TMI = zeros(measurements)
    for measurement in 1:measurements
        for layer in 1:3
            circuit!(state, px, py, pz)
        end
        entropy[measurement, :] = HC_EE(state, L)
        TMI[measurement] = HC_TMI(state, L)
    end

    return sum(entropy, dims=1)./measurements, sum(TMI)./measurements
end

function _chain_trajectory(id, init, L, parameters, mode, checkpoints, verbose)
    
    # Setup
    if mode == :ChainPP
        circuit! = _chain_PP!
    elseif mode == :ChainPQ
        circuit! = _chain_PQ!
    end
    px = parameters[1]
    py = parameters[2]
    pz = parameters[3]
    measurements = 200

    # Check for existing checkpoints, thermalisation
    if checkpoints && isfile("data/checkpoints/$(id).jld2")
        state = _get_checkpoint(id)
        verbose ? (@info "Successfully loaded checkpoint $(id).") : nothing
    else
        (verbose && checkpoints) ? (@info "No checkpoint found for $(id). Starting thermalisation.") : nothing
        state = init
        for layer in 1:3*L
            circuit!(state, px, py, pz)
        end
        verbose ? (@info "Thermalisation complete.") : nothing
        if checkpoints
            _set_checkpoint(id, state)
            verbose ? (@info "Checkpoint $(id) saved.") : nothing
        end
    end

    ### Measurements
    entropy = zeros(measurements, 33)
    TMI = zeros(measurements)
    for measurement in 1:measurements
        for layer in 1:31
            circuit!(state, px, py, pz)
        end
        entropy[measurement, :] = Chain_EE(state)
        TMI[measurement] = Chain_TMI(state)
    end
    
    return sum(entropy, dims=1)./measurements, sum(TMI)./measurements
end