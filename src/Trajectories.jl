import Base: hash, run
abstract type Trajectory end

struct Measurement
    meas_id::Int64
    entropy::Vector{Float64}
    tmi::Float64
    parameters::Vector{Any}
end

"""
    hash(trajectory::Trajectory) -> UInt

    Returns a hash of the trajectory properties.
"""
function hash(trajectory::Trajectory)
    return hash("$(trajectory.name)_$(trajectory.size)_$(trajectory.params)_$(trajectory.index)_$(trajectory.thermalization_steps)_$(trajectory.measurement_steps)_$(trajectory.number_of_measurements)")
end

function _write_checkpoint(state, time, filename, trajectory::Trajectory)
    trajectory.verbosity == :debug ? (@info "Spawned _write_checkpoint subroutine.") : nothing
    jldopen(filename, "w") do file
        println("Writing checkpoint to file $(filename).")
        file["state"] = state
        file["time"] = time
    end
    trajectory.verbosity == :debug ? (@info "Saved checkpoint.") : nothing
end

function _read_checkpoint(filename, trajectory::Trajectory)
    trajectory.verbosity == :debug ? (@info "Sequentially running _read_checkpoint subroutine.") : nothing
    state = jldopen(filename, "r") do file
        file["state"]
    end
    time = jldopen(filename, "r") do file
        file["time"]
    end
    trajectory.verbosity == :debug ? (@info "Loaded checkpoint.") : nothing
    return state, time
end

function checkpoint(state, trajectory::Trajectory, time)
    trajectory.verbosity == :debug ? (@info "Sequentially running checkpoint subroutine.") : nothing

    filename = "data/checkpoints/$(hash(trajectory)).jld2" 
    
    if isfile(filename) # checkpoint file exists
        
        trajectory.verbosity == :debug ? (@info "Found checkpoint file $(filename), checking age...") : nothing
        
        checkpoint_time = jldopen(filename, "r") do file
            file["time"]
        end

        trajectory.verbosity == :debug ? (@info "Checkpoint time is $(checkpoint_time).") : nothing

        if time > checkpoint_time # checkpoint file is older than current time -> overwrite

            trajectory.verbosity == :debug ? (@info "Checkpoint @ $(checkpoint_time), trajectory @ $(time) -> overwrite (async)...") : nothing
            fileio = Threads.@spawn _write_checkpoint(copy(state), copy(time), filename, trajectory)
            
        else # checkpoint file is newer than current time -> load checkpoint

            trajectory.verbosity == :debug ? (@info "Checkpoint @ $(checkpoint_time), trajectory @ $(time) -> load...") : nothing
            
            state, time = _read_checkpoint(filename, trajectory)

            if checkpoint_time != time
                @warn "Checkpoint time $(checkpoint_time) does not match loaded time $(time), something went wrong!"
            else
                trajectory.verbosity == :debug ? (@info "Successfully loaded checkpoint at time $(time).") : nothing
                trajectory.verbosity == :info ? (@info "Successfully loaded checkpoint at time $(time).") : nothing
            end
        end

    else # checkpoint file does not exist

        trajectory.verbosity == :debug ? (@info "No checkpoint file found, creating new checkpoint file...") : nothing
        fileio = Threads.@spawn _write_checkpoint(copy(state), copy(time), filename, trajectory)
        
    end
    return state, time
end

function thermalise!(state, trajectory::Trajectory, time)
    trajectory.verbosity == :debug ? (@info "Thermalising state...") : nothing
    while time < trajectory.thermalization_steps
        if trajectory.checkpoints && time % 100 == 0
            trajectory.verbosity == :debug ? (@info "Calling checkpoint subroutine at time $(time).") : nothing
            state, time = checkpoint(state, trajectory, time)
        end
        circuit!(state, trajectory)
        time += 1
    end
    trajectory.verbosity == :debug ? (@info "Thermalised state at time $(time).") : nothing
end

function observe!(state::QuantumClifford.AbstractStabilizer, trajectory::Trajectory, time)
    for meas_id in 1:trajectory.number_of_measurements
        if trajectory.checkpoints && time % 100 < trajectory.measurement_steps
            trajectory.verbosity == :debug ? (@debug "Calling checkpoint subroutine at time $(time).") : nothing
            state, time = checkpoint(copy(state), trajectory, copy(time))
        end
        for step in 1:trajectory.measurement_steps
            circuit!(state, trajectory)
            time += 1
        end
        trajectory.verbosity == :debug ? (@debug "Spawning (async) measure subroutine at time $(time).") : nothing
        Threads.@spawn measure(copy(state), trajectory, copy(meas_id))
        trajectory.verbosity == :debug ? (@debug "Continue observe routine.") : nothing

    end
end

function measure(state, trajectory::Trajectory, meas_id)
    trajectory.verbosity == :debug ? (@debug "Measuring state (async)...") : nothing
    meas = Measurement(meas_id, entropy(state, trajectory), tmi(state, trajectory), trajectory.params)
    trajectory.verbosity == :debug ? (@debug "Writing measurement $(meas_id) (async) to file...") : nothing
    jldopen("data/measurements/$(hash(trajectory)).jld2", "a+") do file
        file["$(meas_id)"] = meas
    end
    trajectory.verbosity == :debug ? (@debug "Finished measurement $(meas_id) of $(trajectory.number_of_measurements).") : nothing
end

function run(trajectory::Trajectory)
    if !isdir("data/measurements")
        mkdir("data/measurements")
    end
    tick = now()
    time = 0
    trajectory.verbosity == :debug ? (@debug "Starting trajectory $(trajectory.index) of $(trajectory.name).") : nothing
    state = initialise(trajectory)
    trajectory.verbosity == :debug ? (@debug "Initialised state.") : nothing
    thermalise!(state, trajectory, time)
    observe!(state, trajectory, time)
    tock = now()
    trajectory.verbosity == :info ? (@info "Trajectory time: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS")).") : nothing
    
    trajectory.verbosity == :debug ? (@debug "Finished trajectory $(trajectory.index) of $(trajectory.name).") : nothing
    trajectory.verbosity == :debug ? (@info "Time elapsed: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS")).") : nothing
end

