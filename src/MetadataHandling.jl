"""
    archivedTrajectories()

    Returns a vector of all trajectory hashes in the archive.
"""
function archivedTrajectories()
    archivedtrajectories = []
    if isfile("data/archive.jld2")
        jldopen("data/archive.jld2", "r") do file
            for key in keys(file)
                push!(archivedtrajectories, key)
            end
        end
    end
    return parse.(UInt, archivedtrajectories)
end

"""
    boolComplete(traj::Trajectory) :: Bool

    Checks if a trajectory is complete, i.e. if all measurements are present. Automatically averages if all measurements are present but no average is found.
"""
function boolComplete(traj::Trajectory) :: Bool
    if boolArchived(traj)
        traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Trajectory $(hash(traj)) already archived.") : nothing
        return true
    end
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
                writeAverage(traj)
                return true
            end
        end
    else
        traj.verbosity == :debug ? (@info "[TrajectoryComplete?] No file exists for trajectory $(hash(traj)).") : nothing
        return false
    end
end

"""
    boolComplete(sim::Simulation) :: Bool

    Checks if a simulation is complete, i.e. if all linked trajectories are complete.
"""
function boolComplete(sim::Simulation) :: Bool
    archived = archivedTrajectories()
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            if !(hash(sim.ensemble[i][j]) in archived)
                if !boolComplete(sim.ensemble[i][j])
                    return false
                end
            end
        end
    end
    return true
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
    boolArchived(sim::Simulation) :: Bool

    Checks if a simulation is archived, i.e. if all trajectories are archived.
"""
function boolArchived(sim::Simulation) :: Bool
   archived = archivedTrajectories()
    for i in eachindex(sim.ensemble)
          for j in eachindex(sim.ensemble[i])
                if !(hash(sim.ensemble[i][j]) in archived)
                 return false
                end
          end
     end
     return true
end

"""
    archive(traj::Trajectory) :: Nothing

    Archives a trajectory by first checking if it is complete and then archiving it.
"""
function archive(traj::Trajectory) :: Nothing
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
                traj.verbosity == :debug ? (@info "[archive] Trajectory $(hash(traj)) archived.") : nothing
            catch
                @warn "[archive] Unknown error trying to finalize trajectory $(hash(traj))."
            end
        end
    else
         @warn "[archive] Tried to finalize trajectory $(hash(traj)), which has not been completed!."
    end
    return nothing
end

"""
    countTrajectories(sim::Simulation) :: Int

    Counts the number of trajectories in a simulation.
"""
function countTrajectories(sim::Simulation) :: Int
    number_of_trajectories = 0
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            number_of_trajectories += 1
        end
    end
    return number_of_trajectories
end

"""
    archive(sim::Simulation) :: Nothing

    Archives all trajectories in a simulation.
"""
function archive(sim::Simulation) :: Nothing
    if boolArchived(sim)
        @info "Trajectories already archived."
        return nothing
    end
    
    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            archive(sim.ensemble[i][j])
        end
    end

    if boolArchived(sim)
        println("All trajectories archived.")
        println("Do you want to delete the original files? (y/n)")
        answer = readline()
        if answer == "y"
            removeTrajectories(sim)
        end
    else
        println("Something went wrong. Not all trajectories were archived.")
    end
    return nothing
end

"""
    removeTrajectories(sim::Simulation) :: Nothing

    Removes all trajectories in a simulation from the filesystem.
"""
function removeTrajectories(sim::Simulation) :: Nothing
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
    return nothing
end

"""
    writeMetadata(sim::Simulation) :: Nothing

    Writes the metadata of a simulation to the filesystem.
"""
function writeMetadata(sim::Simulation)
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

"""
    loadMetadata()

    Loads all metadata from the filesystem and returns a vector of simulations.
"""
function loadMetadata()
    files = readdir("data/metadata")
    simulations = Vector{Simulation}(undef, length(files))
    for i in eachindex(files)
        file = files[i]
        jldopen("data/metadata/$(file)", "r") do file
           simulations[i] = file["simulation"]
        end
    end
    return simulations
end

"""
    loadMetadata(name::String)

    Loads the metadata of a simulation from the filesystem and returns it. If the simulation name cannot be matched, it loads all metadata.
"""
function loadMetadata(name::String)
   sims = loadMetadata()
    for sim in sims
         if sim.name == name
              return sim
         end
    end
    @info "Simulation $(name) does not exist, loading all metadata..."
    sleep(1)
    return loadMetadata()
end

function metadata()
    sims = loadMetadata()
    println("Found $(length(sims)) simulations.")
    printMetadata(sims)

    num_trajectories = length(allTrajectories())
    println("$(num_trajectories) trajectories.")
    num_adopted = length(adoptedTrajectories())
    num_archived = length(archivedTrajectories())
    if num_adopted == num_trajectories
        println("All trajectories adopted.")
    else
        println("$(num_adopted) adopted trajectories.")
    end
    if num_archived == num_trajectories
        println("All trajectories archived.")
    else
        println("$(num_archived) archived trajectories.")
    end

    num_missing = length(missingTrajectories(sims))
    if  num_missing > 0
        println("$(num_missing) missing trajectories.")
    end

    num_orphaned = length(orphanedTrajectories())
    if num_orphaned > 0
        println("$(num_orphaned) orphaned trajectories.")
    end
end
    
    

function printMetadata(sims::Vector{Simulation}) :: Nothing #custom printing function for simulations vector
    data = String[]
    for sim in sims
        push!(data,"$(sim.name)")
        push!(data,"$(sim.type)")
        push!(data,"$(sim.size)")
        push!(data,"$(length(sim.parameter_set)) × $(sim.num_trajectrories) × $(sim.number_of_measurements)")
        push!(data,"$(boolComplete(sim))")
        push!(data,"$(boolArchived(sim))")
    end
    data = permutedims(reshape(data, 6, length(sims)))
    pretty_table(data; header=["Filename", "Type", "L", "Params. × Traj. × Meas.", "Complete?", "Archived?"])
end

Base.show(io::IO, sims::Vector{Simulation}) = printMetadata(sims) #extend show function with custom pretty printing

archiveTrajectories(string::String) = archiveTrajectories(loadMetadata(string))

archiveTrajectories() = archiveTrajectories.([sim for sim in loadMetadata()])

"""
    allTrajectories()

    Returns a vector of all trajectory hashes in the filesystem.
"""
function allTrajectories() 
    hashes = []
    
    # archived trajectories
    if isfile("data/archive.jld2")
        jldopen("data/archive.jld2", "r") do file
            for key in keys(file)
                push!(hashes, key)
            end
        end
    end

    # unarchived trajectories
    files = readdir("data/measurements")
    for file in files
        push!(hashes, split(file, ".")[1])
    end

    return parse.(UInt, hashes)
end

"""
    orphanedMeasurements()

    Returns a vector of all trajectory hashes in the filesystem that are not in the metadata.
"""
function orphanedTrajectories()
    alltrajectories = allTrajectories()
    adoptedtrajectories = adoptedTrajectories()
    return setdiff(alltrajectories, adoptedtrajectories)
end

"""
    adoptedTrajectories()

    Returns a vector of all trajectory hashes in the filesystem that are in the metadata.
"""
function adoptedTrajectories()
    adoptedtrajectories = []
    archivedtrajectories = archivedTrajectories()
    sims = loadMetadata()
    for sim in sims
        for i in eachindex(sim.ensemble)
            for j in eachindex(sim.ensemble[i])
                if hash(sim.ensemble[i][j]) in archivedtrajectories
                    push!(adoptedtrajectories, hash(sim.ensemble[i][j]))
                elseif isfile("data/measurements/$(hash(sim.ensemble[i][j])).jld2")
                    push!(adoptedtrajectories, hash(sim.ensemble[i][j]))
                end
            end
        end
    end
    return adoptedtrajectories
end

"""
    backup()

    Backs up all metadata and measurements (including the archive) to a folder called "backup" in the data folder.
"""
function backup()
    if !isdir("data/backup")
        mkdir("data/backup")
    end
    if !isdir("data/backup/measurements")
        mkdir("data/backup/measurements")
    end
    if !isdir("data/backup/metadata")
        mkdir("data/backup/metadata")
    end
    try
        if isfile("data/archive.jld2")
            cp("data/archive.jld2", "data/backup/archive.jld2")
        end
        # read measurement files
        files = readdir("data/measurements")
        if !(isempty(files))
            for file in files
                cp("data/measurements/$(file)", "data/backup/measurements/$(file)")
            end
        end
        # read metadata files
        files = readdir("data/metadata")
        if !(isempty(files))
            for file in files
                cp("data/metadata/$(file)", "data/backup/metadata/$(file)")
            end
        end
        @info "Backup successful. Please rename the backup folder to avoid overwriting."
        return true
    catch
        @warn "Backup failed."
        return false
    end
end

"""
    cleanup()

    Cleans up the filesystem by archiving all trajectories in the metadata and deleting all orphaned trajectories.
"""
function cleanup()
    println("Perform backup before cleanup? (y/n)")
    answer = readline()
    backupcomplete = false
    if answer == "y"
        backupcomplete = backup()
    end
    if !backupcomplete
        println("Backup failed. Continue anyway? (y/n)")
        answer = readline()
        if answer != "y"
            println("Aborting.")
            return nothing
        end
    end

    # readout all adopted trajectories from archive
    all = allTrajectories()
    adopted = adoptedTrajectories()
    orphaned = setdiff(all, adopted)
    if isempty(orphaned)
        println("Metadata matches trajectories. Nothing to clean up.")
        return nothing
    end
    adopted_archive = []
    if isfile("data/archive.jld2")
        jldopen("data/archive.jld2", "r") do file
            for key in keys(file)
                push!(adopted_archive, file[key])
            end
        end
    end
    # write new archive
    jldopen("data/archive.jld2", "w") do file
        for i in eachindex(adopted)
            file["$(adopted[i])"] = adopted_archive[i]
        end
    end

    # archive free trajectories
    sims = loadMetadata()
    for sim in sims
        for i in eachindex(sim.ensemble)
            for j in eachindex(sim.ensemble[i])
                if !(hash(sim.ensemble[i][j]) in adopted)
                    archive(sim.ensemble[i][j])
                end
            end
        end
    end

    # delete remaining files
    for hash in orphaned
        if isfile("data/measurements/$(hash).jld2")
            rm("data/measurements/$(hash).jld2")
        end
    end
    println("Cleanup complete.")
    return nothing
end


"""
    missingTrajectories(simulation::Simulation)

    Returns a vector of all trajectories in a simulation that are not complete.
"""
function missingTrajectories(simulation::Simulation)
    missing_traj = []
    if !boolArchived(simulation)
        if !boolComplete(simulation)
            for i in eachindex(simulation.ensemble)
                for j in eachindex(simulation.ensemble[i])
                    if !boolComplete(simulation.ensemble[i][j])
                        push!(missing_traj, simulation.ensemble[i][j])
                    end
                end
            end
        end
    end
    return missing_traj
end

function missingTrajectories(sims::Vector{Simulation})
    missing_traj = []
    for sim in sims
        push!(missing_traj, missingTrajectories(sim))
    end
    return reduce(vcat, missing_traj)
end

missingTrajectories(string::String) = missingTrajectories(loadMetadata(string))

# function mergeArchives(archive1, archive2)
# function backupArchive()

