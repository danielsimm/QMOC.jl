"""
    boolComplete(traj::Trajectory) :: Bool

    Checks if a trajectory is complete, i.e. if all measurements are present. Automatically averages if all measurements are present but no average is found.
"""
function boolComplete(traj::Trajectory) :: Bool
    if isfile("data/measurements/$(hash(traj)).jld2")
        jldopen("data/measurements/$(hash(traj)).jld2", "r") do file
            try # check if average is there
                file["average"]
                traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Found average for trajectory $(hash(traj)).") : nothing
                return true # average found
            catch # no average found
                # check if all measurements are there
                for i in 1:traj.number_of_measurements
                    try 
                        file["$(i)"]
                    catch # measurement missing
                        traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Missing measurement $(i) for trajectory $(hash(traj)).") : nothing
                        return false
                    end
                end
                traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Found all measurements for trajectory $(hash(traj)). Average missing for unknown reasons.") : nothing
                average(traj)
                return true
            end
        end
    else
        traj.verbosity == :debug ? (@info "[TrajectoryComplete?] No file exists for trajectory $(hash(traj)).") : nothing
        return false
    end
end

"""
    boolArchived(traj::Trajectory) :: Bool

    Checks if a trajectory is archived, i.e. if the measurement average is present in the archive.
"""
function boolArchived(traj::Trajectory) :: Bool
    if isfile("data/archive.jld2") # check if archive exists
        jldopen("data/archive.jld2", "r") do file
            try
                file["$(hash(traj))"] # check if trajectory is in archive
                return true
            catch
                return false
            end
        end
    else
        return false
    end
end

"""
    finalize(traj::Trajectory)

    Finalizes a trajectory by first checking if it is complete and then archiving it.
"""
function finalize(traj::Trajectory) 
   if boolArchived(traj)
        traj.verbosity == :debug ? (@info "[finalize] Trajectory $(hash(traj)) already archived.") : nothing
        return nothing
   elseif boolComplete(traj)
        jldopen("data/measurements/$(hash(traj)).jld2", "r") do file
            try
                average = file["average"]
                jldopen("data/archive.jld2", "a+") do file
                    file["$(hash(traj))"] = average
                end
                traj.verbosity == :debug ? (@info "[finalize] Trajectory $(hash(traj)) archived.") : nothing
            catch
                traj.verbosity == :debug ? (@info "[finalize] Trajectory $(hash(traj)) not complete.") : nothing
            end
        end
    else
         @warn "[finalize] Trajectory $(hash(traj)) could not be finalized."
    end
    return nothing
end

"""
    readTrajectory(traj::Trajectory) :: Measurement

    Reads a trajectory from the archive or from the filesystem. Assumes that the trajectory is complete.
"""
function readTrajectory(traj::Trajectory) :: Measurement
    if boolArchived(traj)
        jldopen("data/archive.jld2", "r") do file
            return file["$(hash(traj))"]
        end
    else
        jldopen("data/measurements/$(hash(traj)).jld2", "r") do file
            return file["average"]
        end
    end
end


function commitMetadata(simulation::Simulation)
    if !isdir("data")
        mkdir("data")
    end
    if !isdir("data/metadata")
        mkdir("data/metadata")
    end

    try
        jldopen("data/metadata/$(simulation.name).jld2", "a+") do file
            file["simulation"] = simulation
        end
        push!(SimulationArchive, simulation)
    catch
        @info "Metadata file for simulation $(simulation.name) already exists."
    end

end

function boolArchived(sim::Simulation) :: Bool
    if isfile("data/archive.jld2")
        for i in eachindex(sim.ensemble)
            for j in eachindex(sim.ensemble[i])
               if !boolArchived(sim.ensemble[i][j])
                    return false
                end
            end
        end
    else
        return false
    end
    return true
end

function collectTrajectories(sim::Simulation)
    measurements = []
    hashes = []
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            trajectoryHash = hash(sim.ensemble[i][j])
            if isfile("data/measurements/$(trajectoryHash).jld2")
                push!(hashes, trajectoryHash)
                jldopen("data/measurements/$(trajectoryHash).jld2", "r") do file
                    push!(measurements, file["average"])
                end
            else
                @warn "Trajectory $(trajectoryHash) does not exist."
            end
        end
    end
    return hashes, measurements
end

function archiveTrajectories(sim::Simulation)
    if boolArchived(sim)
        @info "Trajectories already archived."
        return nothing
    end
    
    hashes, measurements = collectTrajectories(sim)

    jldopen("data/archive.jld2", "a+") do file
        for i in eachindex(hashes)
            file["$(hashes[i])"] = measurements[i]
        end
    end

    println("Archived $(length(hashes)) trajectories. Do you want to delete the original files? (y/n)")
    answer = readline()
    if answer == "y"
        removeTrajectories(sim)
    end
end

archiveTrajectories(string::String) = archiveTrajectories(loadSimulation(string))
archiveTrajectories() = archiveTrajectories.([sim for sim in loadSimulationArchive()])


function removeTrajectories(sim::Simulation)
    hashes = []
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            push!(hashes, hash(sim.ensemble[i][j]))
        end
    end
    println("This will delete $(length(hashes)) trajectories from the filesystem. Are you sure? (y/n)")
    answer = readline()
    if answer == "y"
        for hash in hashes
            if isfile("data/measurements/$(hash).jld2")
                rm("data/measurements/$(hash).jld2")
            end
        end
    else
        println("Aborting.")
    end
end


function loadSimulationArchive()
    files = readdir("data/metadata")
    SimulationArchive = []
    for file in files
        jldopen("data/metadata/$(file)", "r") do file
            push!(SimulationArchive, file["simulation"])
        end
    end
    @info "Metadata loaded."
    return SimulationArchive
end


"""
    metadata()

    Prints a table with all metadata after loading into global memory.
"""
function metadata()
    global SimulationArchive = loadSimulationArchive()
    data = String[]
    numberOfSimulations = length(SimulationArchive)
    for simulation in SimulationArchive
        push!(data,"$(simulation.name)")
        push!(data,"$(simulation.type)")
        push!(data,"$(simulation.size)")
        push!(data,"$(simulation.num_trajectrories)")
        push!(data,"$(simulation.number_of_measurements)")
        push!(data,"$(checkMetadataIntegrity(simulation))")
    end
    data = permutedims(reshape(data, 6, numberOfSimulations))
    pretty_table(data; header=["Filename", "Simulation Type", "L", "# Trajectories", "# Measurements", "Metadata OK?"])
end



function loadSimulation(string)
    for simulation in SimulationArchive
        if simulation.name == string
            return simulation
        end
    end
    @info "Simulation $(string) does not exist, loading metadata..."
    sleep(1)
    metadata()
    return nothing
end

function orphanedMeasurements()
    alltrajectories = readdir("data/measurements")
    indexedtrajectories = []
    for sim in loadSimulationArchive()
        for i in eachindex(sim.ensemble)
            for j in eachindex(sim.ensemble[i])
                push!(indexedtrajectories, hash(sim.ensemble[i][j]))
            end
        end
    end
    for indexedtraj in indexedtrajectories
        name = "$(indexedtraj).jld2"
        if name in alltrajectories
            deleteat!(alltrajectories, findfirst(x->x==name, alltrajectories))
        end
    end
    return alltrajectories
end

function missingTrajectories(simulation::Simulation)
    missing_traj = []
    for i in eachindex(simulation.ensemble)
        for j in eachindex(simulation.ensemble[i])
            trajectoryHash = hash(simulation.ensemble[i][j])
            if !isfile("data/measurements/$(trajectoryHash).jld2")
                push!(missing_traj, simulation.ensemble[i][j])
            end
        end
    end
    return missing_traj
end

missingTrajectories(string::String) = missingTrajectories(loadSimulation(string))

