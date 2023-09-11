using PrettyTables

function checkMetadataIntegrity(simulation::Simulation)
    isCorrupted = false
    for i in eachindex(simulation.ensemble)
        for j in eachindex(simulation.ensemble[i])
            trajectoryHash = hash(simulation.ensemble[i][j])
            if !isfile("data/measurements/$(trajectoryHash).jld2")
                @warn "Integrity corrupted: Trajectory $(trajectoryHash) is part of simulation $(simulation.name) from metadata but can not be found at path 'data/measurements/$(trajectoryHash).jld2'."
                isCorrupted = true
            end
        end
    end
    return isCorrupted
end

function checkMetadataIntegrity()
    integrity = true
    if isfile("data/metadata.jld2")
        jldopen("data/metadata.jld2", "r") do file
            for simulationName in keys(file)
                integrity *= !checkMetadataIntegrity(file[simulationName])
            end
        end
    end
    @info "Metadata integrity: $(integrity ? "OK" : "CORRUPTED")"
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

    data = String[]
    jldopen("data/metadata.jld2", "r") do file
        for simulationName in keys(file)
            simulation = file[simulationName]
            push!(data,"$(simulation.name)")
            push!(data,"$(simulation.type)")
            push!(data,"$(simulation.size)")
            push!(data,"$(simulation.num_trajectrories)")
            push!(data,"$(simulation.number_of_measurements)")
        end
    end
    pretty_table(data, ["name", "type", "L", "# trajectories", "# measurements"])
end


