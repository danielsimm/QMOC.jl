
"""
    readTrajectory(traj::Trajectory) :: Measurement

    Reads a trajectory from the archive or from the filesystem. Assumes that the trajectory is complete.
"""
function readTrajectory(traj::Trajectory) :: Measurement
    if boolArchived(traj)
        jldopen("data/archive.jld2", "r") do file
            return file["$(hash(traj))"]
        end
    elseif boolComplete(traj)
        jldopen("data/measurements/$(hash(traj)).jld2", "r") do file
            return file["average"]
        end
    else
        error("Tried to read incomplete trajectory.")
    end
end


function readSimulation(sim::Simulation) :: Array{Measurement}
    N = length(sim.ensemble)
    M = sim.num_trajectrories
    readarray = Array{Measurement}(undef, N, M)
    if boolArchived(sim)
        jldopen("data/archive.jld2", "r") do file
            for i in 1:N
                for j in 1:M
                    readarray[i, j] = file["$(hash(sim.ensemble[i][j]))"]
                end
            end
        end
    else
        for i in 1:N
            for j in 1:M
                readarray[i, j] = readTrajectory(sim.ensemble[i][j])
            end
        end
    end
    return readarray
end

function readSimulation(sim_name::String) :: Array{Measurement}
    return readSimulation(loadMetadata(sim_name))
end

function average(sim::Simulation)
    data = readSimulation(sim)
    N = length(sim.ensemble)
    M = sim.num_trajectrories
    EE = []
    ΔEE = []
    TMI = []
    ΔTMI = []
    for i in 1:N
        push!(EE, mean([data[i, j].entropy for j in 1:M]))
        push!(ΔEE, std([data[i, j].entropy for j in 1:M]))
        push!(TMI, mean([data[i, j].tmi for j in 1:M]))
        push!(ΔTMI, std([data[i, j].tmi for j in 1:M]))
    end
    return EE, ΔEE, TMI, ΔTMI
end


function average(traj_vec::Vector{Trajectory})
    EE_arr = [readTrajectory(traj).entropy for traj in traj_vec]
    TMI_arr = [readTrajectory(traj).tmi for traj in traj_vec]
    EE = mean(permutedims(reduce(hcat, EE_arr)), dims=1)
    ΔEE = std(permutedims(reduce(hcat, EE_arr)), dims=1)
    TMI = mean(TMI_arr)
    ΔTMI = std(TMI_arr)
    return EE, ΔEE, TMI, ΔTMI
end

## TODO: write faster versions for simulations (especially archived ones)
function evaluate(simulation::Simulation, observable::Symbol)
    avg = average(simulation)
    if observable == :I3
        return simulation.parameter_set, avg[3]
    elseif observable == :ΔI3
        return simulation.parameter_set, avg[4]
    elseif observable == :SvN
        return simulation.parameter_set, avg[1], subsystem_labels(simulation.ensemble[1][1])
    elseif observable == :ΔSvN
        return simulation.parameter_set, avg[2], subsystem_labels(simulation.ensemble[1][1])
    else
        error("Observable $(observable) not implemented, choose from [:I3, :ΔI3, :SvN, :ΔSvN].")
    end
end

evaluate(sim_name::String, Observable::Symbol) = evaluate(loadMetadata(sim_name), Observable)
evaluate(simulation::Simulation) = evaluate(simulation, :SvN)
evaluate(sim_name::String) = evaluate(sim_name, :SvN)