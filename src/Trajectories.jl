import Base: hash, run
abstract type Trajectory end

# testtraj = KitaevTrajectory(4, 10, "test", [1//2, 1//2, 1//2], true, :high, 4, 100, 100, 100)

struct Measurement
    meas_id::Int64
    entropy::Vector{Float64}
    tmi::Float64
    parameters::Vector{Real}
end

"""
    hash(trajectory::Trajectory) -> UInt

    Returns a hash of the trajectory properties.
"""
function hash(trajectory::Trajectory)
    return hash("$(trajectory.name)_$(trajectory.size)_$(trajectory.params)_$(trajectory.index)_$(trajectory.thermalization_steps)_$(trajectory.measurement_steps)_$(trajectory.number_of_measurements)")
end

function checkpoint(state, trajectory::Trajectory, time)
    filename = "data/checkpoints/$(hash(trajectory)).jld2" 
    if isfile(filename) # checkpoint file exists
        trajectory.verbosity == :high ? (@info "Found checkpoint file $(filename), checking age...") : nothing
        checkpoint_time = jldopen(filename, "r") do file
            file["time"]
        end
        if time > checkpoint_time # checkpoint file is older than current time -> overwrite
            trajectory.verbosity == :high ? (@info "Checkpoint @ $(checkpoint_time), trajectory @ $(time) -> overwrite...") : nothing
            fileio = Threads.@spawn begin
                jldopen(filename, "w") do file # overwrite checkpoint file
                    file["state"] = state
                    file["time"] = time
                end
                trajectory.verbosity == :high ? (@info "Done!") : nothing
            end
        else # checkpoint file is newer than current time -> load checkpoint
            trajectory.verbosity == :high ? (@info "Checkpoint @ $(checkpoint_time), trajectory @ $(time) -> load...") : nothing
            state = jldopen(filename, "r") do file
                file["state"]
            end
            time = jldopen(filename, "r") do file
                file["time"]
            end
            if checkpoint_time != time
                @warn "Checkpoint time $(checkpoint_time) does not match loaded time $(time), something went wrong!"
            else
                trajectory.verbosity == :high ? (@info "Successfully loaded checkpoint at time $(time).") : nothing
                trajectory.verbosity == :low ? (@info "Successfully loaded checkpoint at time $(time).") : nothing
            end
        end
    else # checkpoint file does not exist
        trajectory.verbosity == :high ? (@info "No checkpoint file found, creating new checkpoint file...") : nothing
        fileio = Threads.@spawn begin
            jldopen(filename, "w") do file # create checkpoint file
                file["state"] = state
                file["time"] = time
            end
            trajectory.verbosity == :high ? (@info "Done!") : nothing
        end
    end
    return state, time
    wait(fileio)
end

function thermalize!(state, trajectory::Trajectory, time)
    while time < trajectory.thermalization_steps
        if trajectory.checkpoints && time % 100 == 0
            state, time = checkpoint(copy(state), trajectory, copy(time))
        end
        circuit!(state, trajectory)
        time += 1
    end
end

function observe!(state::QuantumClifford.AbstractStabilizer, trajectory::Trajectory, time)
    for meas_id in 1:trajectory.number_of_measurements
        if trajectory.checkpoints && time % 100 == 0
            state, time = checkpoint(copy(state), trajectory, copy(time))
        end
        for step in 1:trajectory.measurement_steps
            circuit!(state, trajectory)
            time += 1
        end
        Threads.@spawn measure(copy(state), trajectory, copy(meas_id))
    end
    return measurements
end

function measure(state, trajectory::Trajectory, meas_id)
    meas = Measurement(meas_id, entropy(state, trajectory), tmi(state, trajectory), trajectory.params)
    jldopen("data/measurements/$(hash(trajectory)).jld2", "a+") do file
        file["$(meas_id)"] = meas
    end
end

function run(trajectory::Trajectory)
    time = 0
    state = initialise(trajectory)
    thermalize!(state, trajectory, time)
    observe!(state, trajectory, time)
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