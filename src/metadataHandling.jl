function checkMetadataIntegrity(simulation::Simulation)
    isCorrupted = false
    for i in eachindex(simulation.ensemble)
        for j in eachindex(simulation.ensemble[i])
            trajectoryHash = hash(simulation.ensemble[i][j])
            if !isfile("data/measurements/$(trajectoryHash).jld2")
                @warn "Integrity corrupted: Expecting trajectory $(trajectoryHash) as part of simulation $(simulation.name) from metadata but can not be found at path 'data/measurements/$(trajectoryHash).jld2'."
                isCorrupted = true
            end
        end
    end
    return isCorrupted
end

function checkMetadataIntegrity()
    integrity = true
    files = readdir("data/metadata")
    if !(isempty(files))
        for file in files
            jldopen("data/metadata/$(file)", "r") do file
                integrity *= !checkMetadataIntegrity(file["simulation"])
            end
        end
    end
    @info "Metadata integrity: $(integrity ? "OK" : "CORRUPTED")"
end

function commitMetadata(simulation::Simulation)
    if !isdir("data")
        mkdir("data")
    end
    if !isdir("data/metadata")
        mkdir("data/metadata")
    end

    jldopen("data/metadata/$(simulation.name).jld2", "a+") do file
        file["simulation"] = simulation
    end

    push!(SimulationArchive, simulation)
end

function printMetadata()
    # files = readdir("data/metadata")
    # if isempty(files)
    #     @warn "No metadata file found at path 'data/metadata.jld2'."
    #     return
    # end
    # data = String[]
    # numberOfSimulations = length(files)
    # for file in files
    #     jldopen(file, "r") do file
    #         simulation = file[simulation]
    #         push!(data,"$(simulation.name)")
    #         push!(data,"$(simulation.type)")
    #         push!(data,"$(simulation.size)")
    #         push!(data,"$(simulation.num_trajectrories)")
    #         push!(data,"$(simulation.number_of_measurements)")
    #     end
    # end
    data = String[]
    numberOfSimulations = length(SimulationArchive)
    for simulation in SimulationArchive
        push!(data,"$(simulation.name)")
        push!(data,"$(simulation.type)")
        push!(data,"$(simulation.size)")
        push!(data,"$(simulation.num_trajectrories)")
        push!(data,"$(simulation.number_of_measurements)")
    end
    data = permutedims(reshape(data, 5, numberOfSimulations))
    pretty_table(data; header=["Filename", "Simulation Type", "L", "Number of Trajectories", "Number of Measurements"])
    
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

function loadSimulation(string)
    for simulation in SimulationArchive
        if simulation.name == string
            return simulation
        end
    end
    @info "Simulation $(string) does not exist, printing metadata..."
    sleep(1)
    printMetadata()
    return nothing
end


