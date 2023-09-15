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
    if isfile("data/metadata.jld2")
        jldopen("data/metadata.jld2", "r") do file
            for simulationName in keys(file)
                METADATA_INTEGRITY *= !checkMetadataIntegrity(file[simulationName])
            end
        end
    end
    @info "Metadata integrity: $(METADATA_INTEGRITY ? "OK" : "CORRUPTED")"
end

function commitMetadata(simulation::Simulation)
    if !isdir("data")
        mkdir("data")
    end

    jldopen("data/metadata.jld2", "a+") do file
        file[simulation.name] = simulation
    end
end

function printMetadata()
    
    if !isfile("data/metadata.jld2")
        @warn "No metadata file found at path 'data/metadata.jld2'."
        return
    end
    numberOfSimulations = 0
    data = String[]
    jldopen("data/metadata.jld2", "r") do file
        for simulationName in keys(file)
            numberOfSimulations += 1
            simulation = file[simulationName]
            push!(data,"$(simulation.name)")
            push!(data,"$(simulation.type)")
            push!(data,"$(simulation.size)")
            push!(data,"$(simulation.num_trajectrories)")
            push!(data,"$(simulation.number_of_measurements)")
        end
    end
    reshape(data, 5, numberOfSimulations)
    pretty_table(data, ["name", "type", "L", "# trajectories", "# measurements"])
end

function average_observable(trajectory)
end

function loadSimulationArchive() :: Vector{Simulation}
    SimulationArchive = []
    # load all available Simulations into the global SimulationArchive object
    return SimulationArchive
end

