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
    jldopen(filename, "w") do file
        trajectory.verbosity == :debug ? (@info "Writing checkpoint to file $(filename).") : nothing
        file["state"] = state
        file["time"] = time
    end
    trajectory.verbosity == :debug ? (@info "Closed file $(filename).") : nothing
end

function _read_checkpoint(filename, trajectory::Trajectory)
    trajectory.verbosity == :debug ? (@info "Running _read_checkpoint subroutine.") : nothing
    state = jldopen(filename, "r") do file
        file["state"]
    end
    time = jldopen(filename, "r") do file
        file["time"]
    end
    trajectory.verbosity == :debug ? (@info "Loaded checkpoint.") : nothing
    return state, time
end

function checkpoint(state, trajectory::Trajectory, time, existing_measurements)
    filename = "data/checkpoints/$(hash(trajectory)).jld2" 

    if time == 0 # begining of trajectory, check for previous data
        trajectory.verbosity == :debug ? (@info "Beginning of trajectory, checking for previous data...") : nothing
        if isfile(filename)
            trajectory.verbosity == :debug ? (@info "Found previous data, loading...") : nothing

            state, time = _read_checkpoint(filename, trajectory)
            
            for i in 1:trajectory.number_of_measurements
                try 
                    jldopen(filename, "r") do file
                        file["$(i)"]
                    end
                    existing_measurements += 1
                catch
                    nothing
                end
            end

            trajectory.verbosity == :debug ? (@info "--- Checkpoint info --- \n Time: $(time) \n Thermalised: $(time>trajectory.thermalization_steps) \n Existing measurements: $(existing_measurements) of $(trajectory.number_of_measurements)") : nothing
        else
            trajectory.verbosity == :debug ? (@info "No previous data found, starting from scratch.") : nothing
        end

    else
        _write_checkpoint(copy(state), copy(time), filename, trajectory)
        trajectory.verbosity == :debug ? (@info "--- Checkpoint info --- \n Time: $(time) \n Thermalised: $(time>trajectory.thermalization_steps) \n Existing measurements: $(existing_measurements) of $(trajectory.number_of_measurements)") : nothing
    end
    trajectory.verbosity == :debug ? (@info "Finished checkpoint subroutine.") : nothing
    return state, time, existing_measurements
    
end

function thermalise(state, trajectory::Trajectory, time, existing_measurements)
    trajectory.verbosity == :debug ? (@info "Thermalising state...") : nothing
    while time < trajectory.thermalization_steps
        if trajectory.checkpoints && time % 10 == 0
            trajectory.verbosity == :debug ? (@info "Calling checkpoint subroutine at time $(time).") : nothing
            state, time, existing_measurements = checkpoint(state, trajectory, time, existing_measurements)
        end
        circuit!(state, trajectory)
        time += 1
    end
    trajectory.verbosity == :debug ? (@info "Thermalised state at time $(time).") : nothing
    return state, time, existing_measurements
end

function observe(state::QuantumClifford.AbstractStabilizer, trajectory::Trajectory, time, existing_measurements)
    for meas_id in existing_measurements+1:trajectory.number_of_measurements
        if trajectory.checkpoints && meas_id % 10 == 0
            trajectory.verbosity == :debug ? (@info "Calling checkpoint subroutine at time $(time).") : nothing
            state, time, existing_measurements = checkpoint(state, trajectory, time, existing_measurements)
        end
        for step in 1:trajectory.measurement_steps
            circuit!(state, trajectory)
            time += 1
        end
        trajectory.verbosity == :debug ? (@info "Spawning (async) measure subroutine at time $(time).") : nothing
        Threads.@spawn measure(copy(state), trajectory, copy(meas_id))
        existing_measurements += 1
        trajectory.verbosity == :debug ? (@info "Continue observe routine.") : nothing
    end
    return state, time
end

function measure(state, trajectory::Trajectory, meas_id)
    filename = "data/measurements/$(hash(trajectory)).jld2" 
    trajectory.verbosity == :debug ? (@info "Measurement $(meas_id) subroutine") : nothing
    meas = Measurement(meas_id, entropy(state, trajectory), tmi(state, trajectory), trajectory.params)
    jldopen(filename, "a+") do file
        println("Writing measurement $(meas_id) to file $(filename).")
        file["$(meas_id)"] = meas
    end
    trajectory.verbosity == :debug ? (@info "Finished measurement $(meas_id) of $(trajectory.number_of_measurements).") : nothing
end

function run(trajectory::Trajectory)
    if !isdir("data")
        mkdir("data")
    end
    if !isdir("data/measurements")
        mkdir("data/measurements")
    end
    if !isdir("data/checkpoints")
        mkdir("data/checkpoints")
    end

    trajectory.verbosity == :debug ? (@info "Starting trajectory $(trajectory.index) of $(trajectory.name).") : nothing
    tick = now()
    time = 0
    state = initialise(trajectory)
    trajectory.verbosity == :debug ? (@info "Initialised state.") : nothing

    state, time, existing_measurements = checkpoint(state, trajectory, time, 0)
    
    state, time = thermalise(state, trajectory, time, existing_measurements)
    state, time = observe(state, trajectory, time, existing_measurements)
    tock = now()
    in(trajectory.verbosity, [:info, :debug]) ? (@info "Trajectory time: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS")).") : nothing
    
    trajectory.verbosity == :debug ? (@info "Finished trajectory $(trajectory.index) of $(trajectory.name).") : nothing
end

