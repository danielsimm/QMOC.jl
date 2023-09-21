import Base: hash, run
abstract type Trajectory end

struct Measurement
    entropy::Vector{Float64}
    tmi::Float64
end


"""
    hash(trajectory::Trajectory) -> UInt

    Returns a hash of the trajectory properties.
"""
function hash(trajectory::Trajectory)
    return hash("$(typeof(trajectory))_$(trajectory.size)_$(trajectory.params)_$(trajectory.index)_$(trajectory.thermalization_steps)_$(trajectory.measurement_steps)_$(trajectory.number_of_measurements)")
end

function write_checkpoint(state, time, trajectory::Trajectory)
    filename = "data/checkpoints/$(hash(trajectory)).jld2" 
    jldopen(filename, "w") do file
        trajectory.verbosity == :debug ? (@info "[writecheckpoint] Writing checkpoint @ $(time) to file $(filename).") : nothing
        file["state"] = state
        file["time"] = time
    end
    trajectory.verbosity == :debug ? (@info "[writecheckpoint] Closed file $(filename).") : nothing
end

function read_checkpoint(trajectory::Trajectory)
    filename = "data/checkpoints/$(hash(trajectory)).jld2" 
    state = jldopen(filename, "r") do file
        file["state"]
    end
    time = jldopen(filename, "r") do file
        file["time"]
    end
    trajectory.verbosity == :debug ? (@info "[readcheckpoint] Loaded checkpoint.") : nothing
    return state, time
end

function skippable(trajectory::Trajectory)
    filename2 = "data/measurements/$(hash(trajectory)).jld2"
    existing_measurements = 0
    try
        jldopen(filename2, "r") do file
            file["average"]
        end
        return true
    catch
        if isfile(filename2)
            for i in 1:trajectory.number_of_measurements
                try 
                    jldopen(filename2, "r") do file
                        file["$(i)"]
                    end
                    existing_measurements += 1
                catch
                    nothing
                end
            end
        end
        if existing_measurements == trajectory.number_of_measurements
            return true
            trajectory.verbosity == :debug ? (@info "[skippable] Trajectory $(trajectory.index) of $(trajectory.name) already done -> skip.") : nothing
        else
            return false
        end
    end
end

function hot_start(trajectory::Trajectory)
    filename = "data/checkpoints/$(hash(trajectory)).jld2" 
    filename2 = "data/measurements/$(hash(trajectory)).jld2"
    existing_measurements = 0
    time = 0
    trajectory.verbosity == :debug ? (@info "[hotstart] Attempting hot start, checking for previous data...") : nothing
    if isfile(filename)
        trajectory.verbosity == :debug ? (@info "[hotstart] Found previous data, loading...") : nothing

        state, time = read_checkpoint(trajectory)
        if isfile(filename2)
            trajectory.verbosity == :debug ? (@info "[hotstart] Found previous measurements, counting...") : nothing
            for i in 1:trajectory.number_of_measurements
                try 
                    jldopen(filename2, "r") do file
                        file["$(i)"]
                    end
                    existing_measurements += 1
                catch
                    nothing
                end
            end
        end
        trajectory.verbosity == :debug ? (@info "--- [hotstart] Checkpoint info --- \n Time: $(time) \n Thermalised: $(time>trajectory.thermalization_steps) \n Existing measurements: $(existing_measurements) of $(trajectory.number_of_measurements)") : nothing
    else
        trajectory.verbosity == :debug ? (@info "[hotstart] No previous data found, initialising state...") : nothing
        state = initialise(trajectory)
        trajectory.verbosity == :debug ? (@info "[hotstart] Done.") : nothing
    end
    return state, time, existing_measurements
end

function thermalise(state, trajectory::Trajectory, time, operators::Array{PauliOperator})
    trajectory.verbosity == :debug ? (@info "[thermalise] Thermalising state...") : nothing
    while time < trajectory.thermalization_steps
        circuit!(state, trajectory, operators)
        time += 1
        if trajectory.checkpoints && time % 100 == 0
            trajectory.verbosity == :debug ? (@info "[thermalise] Calling checkpoint subroutine at time $(time).") : nothing
            write_checkpoint(state, time, trajectory)
        end
    end
    trajectory.verbosity == :debug ? (@info "[thremalise] Thermalised state at time $(time).") : nothing
    return state, time
end

function measurement_evolve!(state, trajectory::Trajectory, operators::Array{PauliOperator})
    tick = now()

    for step in 1:trajectory.measurement_steps
        circuit!(state, trajectory, operators)
    end

    tock = now()
    trajectory.verbosity == :debug ? (@info "[measurement_evolve] Finished evolution in $(tock - tick).") : nothing
    return state
end

function measure(state, trajectory::Trajectory, meas_id, time)
    tick = now()

    filename = "data/measurements/$(hash(trajectory)).jld2" 

    trajectory.verbosity == :debug ? (@info "[measure] Measuring data point $(meas_id) of $(trajectory.number_of_measurements...)") : nothing
    meas = Measurement(entropy(state, trajectory), tmi(state, trajectory))

    jldopen(filename, "a+") do file
        trajectory.verbosity == :debug ? (@info "[measure] Writing measurement $(meas_id) to file $(filename).") : nothing
        file["$(meas_id)"] = meas
    end

    # write checkpoint after each measurement
    if trajectory.checkpoints
        trajectory.verbosity == :debug ? (@info "[measure] Calling checkpoint subroutine at time $(time).") : nothing
        write_checkpoint(state, time, trajectory)
    end
    tock = now()

    trajectory.verbosity == :debug ? (@info "[measure] Finished measurement $(meas_id) of $(trajectory.number_of_measurements) in $(tock - tick).") : nothing
end

function observe(state::QuantumClifford.AbstractStabilizer, trajectory::Trajectory, time, existing_measurements, operators)
    for meas_id in existing_measurements+1:trajectory.number_of_measurements-1
        # independently spawn measurement and evolution
        
        @sync begin
            trajectory.verbosity == :debug ? (@info "[observe] Asynchronous measurement $(meas_id) at time $(time).") : nothing
            Threads.@spawn measure(copy(state), trajectory, copy(meas_id), copy(time))
            Threads.@spawn measurement_evolve!(state, trajectory, operators)
            trajectory.verbosity == :debug ? (@info "[observe] Waiting for sync...") : nothing
        end
        trajectory.verbosity == :debug ? (@info "[observe] Synced.") : nothing
        time += trajectory.measurement_steps
        
    end
    if existing_measurements < trajectory.number_of_measurements
        measure(state, trajectory, trajectory.number_of_measurements, time)
    end
    return state, time
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
    if skippable(trajectory)
        return nothing
    end
    trajectory.verbosity == :debug ? (@info "[run] Starting trajectory $(trajectory.index) of $(trajectory.name).") : nothing
    tick = now()
    time = 0

    state, time, existing_measurements = hot_start(trajectory)

    operators = get_operators(trajectory)
    
    state, time = thermalise(state, trajectory, time, operators)
    state, time = observe(state, trajectory, time, existing_measurements, operators)
    
    tock = now()

    trajectory.verbosity == :info ? (@info "[run] Trajectory time: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS")).") : nothing
    trajectory.verbosity == :debug ? (@info "[run] Trajectory time: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS")), trajectory $(trajectory.index) of $(trajectory.name).") : nothing
end

