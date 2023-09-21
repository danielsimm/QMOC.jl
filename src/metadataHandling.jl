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

# function checkMetadataIntegrity()
#     files = readdir("data/metadata")
#     corruptedSimulations = []
#     if !(isempty(files))
#         for file in files
#             jldopen("data/metadata/$(file)", "r") do file
#                 if checkMetadataIntegrity(file["simulation"])
#                     push!(corruptedSimulations, file["simulation"])
#                 end
#             end
#         end
#     end
#     integrity = isempty(corruptedSimulations)
#     if integrity
#         @info "Metadata integrity: OK"
#     else
#         @warn "Metadata integrity: CORRUPTED"
#         @info "Corrupted simulations:"
#         for simulation in corruptedSimulations
#             @info "  $(simulation.name)"
#         end
#         @info "Check for missing trajectories with missingTrajectories(simulation) or missingTrajectories(simulation.name)."
#     end
# end

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

function printMetadata()
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

function removeOrphanedMeasurements()
    # test if all measurements are still present in metadata, remove if not accessible via simulation
    # DO NOT export
    # orphaned measurements get adopted automatically when running a simulation
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

