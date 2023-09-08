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

function thermalise!(state, trajectory::Trajectory, time)
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
    thermalise!(state, trajectory, time)
    observe!(state, trajectory, time)
end

export run, Trajectory, Measurement