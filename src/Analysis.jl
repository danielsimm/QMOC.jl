
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

function average(traj_vec::Vector{Trajectory})
    EE_arr = [average(traj).entropy for traj in traj_vec]
    TMI_arr = [average(traj).tmi for traj in traj_vec]
    EE = mean(permutedims(reduce(hcat, EE_arr)), dims=1)
    ΔEE = std(permutedims(reduce(hcat, EE_arr)), dims=1)
    TMI = mean(TMI_arr)
    ΔTMI = std(TMI_arr)
    return EE, ΔEE, TMI, ΔTMI
end

function evaluate(simulation::Simulation, observable::Symbol)
    if observable == :I3
        return simulation.parameter_set, [average(vec)[3] for vec in simulation.ensemble]
    elseif observable == :ΔI3
        return simulation.parameter_set, [average(vec)[4] for vec in simulation.ensemble]
    elseif observable == :SvN
        return simulation.parameter_set, [average(vec)[1] for vec in simulation.ensemble], subsystem_labels(simulation.ensemble[1][1])
    elseif observable == :ΔSvN
        return simulation.parameter_set, [average(vec)[2] for vec in simulation.ensemble], subsystem_labels(simulation.ensemble[1][1])
    else
        error("Observable $(observable) not implemented, choose from [:I3, :ΔI3, :SvN, :ΔSvN].")
    end
end

evaluate(sim_name::String, Observable::Symbol) = evaluate(loadSimulation(sim_name), Observable)
evaluate(simulation::Simulation) = evaluate(simulation, :SvN)
evaluate(sim_name::String) = evaluate(sim_name, :SvN)