"""
    archivedTrajectories()

    Returns a vector of all trajectory hashes in the archive.
"""
function archivedTrajectories() ::Vector{UInt}
    if isfile("data/archive.jld2")
        ks = jldopen("data/archive.jld2", "r") do file
            keys(file)
        end
        return parse.(UInt, ks)
    else
        return Vector{UInt}(undef, 0)
    end
end

function hashes(sim::Simulation)
    return [hash(sim.ensemble[i][j]) for i in eachindex(sim.ensemble) for j in eachindex(sim.ensemble[i])]
end

"""
    boolComplete(traj::Trajectory) :: Bool

    Checks if a trajectory is complete, i.e. if all measurements are present. Automatically averages if all measurements are present but no average is found.
"""
function boolComplete(traj::Trajectory)::Bool
    if isfile("data/measurements/$(hash(traj)).jld2")
        bool = false
        jldopen("data/measurements/$(hash(traj)).jld2", "r") do file
            if haskey(file, "average") # check if average is there
                # traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Found average for trajectory $(hash(traj)).") : nothing
                bool = true # average found
            else # no average found
                # check if all measurements are there
                bool = true
                for i in 1:traj.number_of_measurements
                    if !(haskey(file, "$(i)")) # measurement missing
                        # traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Missing measurement $(i) for trajectory $(hash(traj)).") : nothing
                        bool = false
                    end
                end
                # traj.verbosity == :debug ? (@info "[TrajectoryComplete?] Found all measurements for trajectory $(hash(traj)). Average missing for unknown reasons.") : nothing
            end
        end
        if bool
            writeAverage(traj)
        end
        return bool
    else
        # traj.verbosity == :debug ? (@info "[TrajectoryComplete?] No file exists for trajectory $(hash(traj)).") : nothing
        return false
    end
end


"""
    boolComplete(sim::Simulation) :: Bool

    Checks if a simulation is complete, i.e. if all linked trajectories are complete.
"""
function boolComplete(sim::Simulation)::Bool
    archived = archivedTrajectories()
    sim_hashes = hashes(sim)
    diff = setdiff(sim_hashes, archived) # non-archived trajectories
    completed = 0

    for i in eachindex(sim.ensemble)
        for j in eachindex(sim.ensemble[i])
            if (hash(sim.ensemble[i][j]) in diff)
                if boolComplete(sim.ensemble[i][j]) 
                    completed += 1
                end
            end
        end
    end

    return (completed == length(diff))
end


"""
    boolArchived(traj::Trajectory) :: Bool

    Checks if a trajectory is archived, i.e. if the measurement average is present in the archive.
"""
function boolArchived(traj::Trajectory)::Bool
    archived = archivedTrajectories()
    traj_hash = hash(traj)
    diff = setdiff(traj_hash, archived)
    return isempty(diff)
end

"""
    boolArchived(sim::Simulation) :: Bool

    Checks if a simulation is archived, i.e. if all trajectories are archived.
"""
function boolArchived(sim::Simulation)::Bool
    archived = archivedTrajectories()
    sim_hashes = hashes(sim)
    diff = setdiff(sim_hashes, archived)
    return isempty(diff)
end

"""
    archive(traj::Trajectory) :: Nothing

    Archives a trajectory by first checking if it is complete and then archiving it.
"""
function archive(traj::Trajectory)::Nothing
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
function countTrajectories(sim::Simulation)::Int
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
function archive(sim::Simulation)::Nothing
    if boolArchived(sim)
        @info "Trajectories already archived."
        return nothing
    end
    
    # check if all trajectories are complete
    if !boolComplete(sim)
        @warn "Not all trajectories are complete. Aborting."
        return nothing
    end

    # collect non-archived hashes
    archived_hashes = archivedTrajectories()
    simulation_hashes = hashes(sim)
    unarchived_hashes = setdiff(simulation_hashes, archived_hashes)
    measurements = Vector{Measurement}(undef, length(unarchived_hashes))
    
    Threads.@threads for i in eachindex(unarchived_hashes)
        jldopen("data/measurements/$(unarchived_hashes[i]).jld2", "r") do file
            measurements[i] = file["average"]
        end
    end
    
    # write to archive
    jldopen("data/archive.jld2", "a+") do file
        for i in eachindex(unarchived_hashes)
            file["$(unarchived_hashes[i])"] = measurements[i]
        end
    end

    if boolArchived(sim)
        println("All trajectories archived.")
        # println("Do you want to delete the original files? (y/n)")
        # answer = readline()
        # if answer == "y"
        #     removeTrajectories(sim)
        # end
    else
        println("Something went wrong. Not all trajectories were archived.")
    end
    return nothing
end

archive(string::String) = archive(loadMetadata(string))

archive() = archive.([sim for sim in loadMetadata()])

"""
    removeTrajectories(sim::Simulation) :: Nothing

    Removes all trajectories in a simulation from the filesystem.
"""
function removeTrajectories(sim::Simulation)::Nothing
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
        jldopen("data/metadata/$(sim.name).jld2", "a+") do file
            file["simulation"] = sim
        end
        @info "Metadata file for simulation $(sim.name) written to disk."
    catch
        @info "Metadata file for simulation $(sim.name) already exists."
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
    if num_missing > 0
        println("$(num_missing) missing trajectories.")
    end

    num_orphaned = length(orphanedTrajectories())
    if num_orphaned > 0
        println("$(num_orphaned) orphaned trajectories.")
    end
end

function printMetadata(sims::Vector{Simulation})::Nothing #custom printing function for simulations vector
    data = String[]
    for sim in sims
        push!(data, "$(sim.name)")
        push!(data, "$(sim.type)")
        push!(data, "$(sim.size)")
        push!(data, "$(length(sim.parameter_set)) × $(sim.num_trajectrories) × $(sim.number_of_measurements)")
        push!(data, "$(boolComplete(sim))")
        push!(data, "$(boolArchived(sim))")
    end
    data = permutedims(reshape(data, 6, length(sims)))
    pretty_table(data; header=["Filename", "Type", "L", "Params. × Traj. × Meas.", "Complete?", "Archived?"])
end

Base.show(io::IO, sims::Vector{Simulation}) = printMetadata(sims) #extend show function with custom pretty printing

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
function cleanup(;interactive::Bool=true)::Nothing
    if interactive
        println("Perform backup before cleanup? (y/n)")
        answer = readline()
        backupcomplete = false
        if answer == "y"
            backupcomplete = backup()
            if !backupcomplete
                println("Backup failed. Continue anyway? (y/n)")
                answer = readline()
                if answer != "y"
                    println("Aborting.")
                    return nothing
                end
            end
        end
    end

    # readout all adopted trajectories from archive
    simulations = loadMetadata()
    all_hashes = allTrajectories()
    adopted_hashes = adoptedTrajectories()
    orphaned_hashes = setdiff(all_hashes, adopted_hashes)
    
    Archive = load_archive()
    # add adopted, non-archived trajectories to archive
    for sim in simulations
        for i in eachindex(sim.ensemble)
            for j in eachindex(sim.ensemble[i])
                if !(hash(sim.ensemble[i][j]) in archived_hashes)
                    if !(hash(sim.ensemble[i][j]) in orphaned_hashes)
                        jldopen("data/measurements/$(hash(sim.ensemble[i][j])).jld2", "r") do file
                            push!(archived_measurements, file["average"])
                        end
                        push!(archived_hashes, hash(sim.ensemble[i][j]))
                    end
                end
            end
        end
    end
    # write new archive
    println("Writing new archive with $(length(archived_hashes)) trajectories... ")
    jldopen("data/archive.jld2", "w") do file
        for i in eachindex(archived_hashes)
            file["$(archived_hashes[i])"] = archived_measurements[i]
        end
    end
    print("Done!")

    # remove archived trajectories
    println("Deleting $(length(archived_hashes)) trajectory files... ")
    for hash in archived_hashes
        if isfile("data/measurements/$(hash).jld2")
            rm("data/measurements/$(hash).jld2")
        end
    end
    print("Done!")

    # delete remaining files
    if interactive
        println("Do you want to delete all orphaned trajectories? (y/n)")
        answer = readline()
        if answer == "y"
            println("Deleting $(length(orphaned_hashes)) orphaned trajectories... ")
            for hash in orphaned_hashes
                if isfile("data/measurements/$(hash).jld2")
                    rm("data/measurements/$(hash).jld2")
                end
            end
            print("Done!")
        end
    end
    

    # delete checkpoints
    files = readdir("data/checkpoints")
    if !(isempty(files))
        for file in files
            rm("data/checkpoints/$(file)")
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


function mergeArchives(archive1, archive2)
    println("Create backup? (y/n)")
    answer = readline()
    if answer == "y"
        backupArchive()
    end
    arch1 = load_archive(archive1)
    arch2 = load_archive(archive2)
    merged = merge(arch1, arch2)
    rm(archive1)
    rm(archive2)
    save(archive1, merged)
end

mergeArchive(path) = mergeArchives("data/archive.jld2", path)

function backupArchive()
    if !isdir("data/backup")
        mkdir("data/backup")
    end
    if !isdir("data/backup/archive")
        mkdir("data/backup/archive")
    end
    try
        if isfile("data/archive.jld2")
            cp("data/archive.jld2", "data/backup/archive/archive_$(today()).jld2")
        end
        @info "Archive backup successful."
    catch
        @warn "Archive backup failed."
    end
end

function load_archive(path)
    return load(path)
end

function load_archive()
    return load_archive("data/archive.jld2")
end
