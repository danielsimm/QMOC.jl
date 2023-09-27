function checkMetadataIntegrity(simulation::Simulation)
    isCorrupted = false
    for i in eachindex(simulation.ensemble)
        for j in eachindex(simulation.ensemble[i])
            trajectoryHash = hash(simulation.ensemble[i][j])
            if !isfile("data/measurements/$(trajectoryHash).jld2")
                #@warn "Integrity corrupted: Expecting trajectory $(trajectoryHash) as part of simulation $(simulation.name) from metadata but can not be found at path 'data/measurements/$(trajectoryHash).jld2'."
                isCorrupted = true
            end
        end
    end
    return isCorrupted
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

function archiveTrajectories(sim::Simulation)
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
            end
        end
    end
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
        push!(data,"$(!checkMetadataIntegrity(simulation))")
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

