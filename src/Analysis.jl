function average(trajectory::Trajectory)
    filename = "data/measurements/$(hash(trajectory)).jld2"
    measurements = Vector{Measurement}(undef, trajectory.number_of_measurements)
    try
        jldopen(filename, "r") do file
            return file["average"]
        end
    catch
        jldopen(filename, "r") do file
            for t in 1:trajectory.number_of_measurements
                measurements[t] = file["$(t)"]
            end
        end
        EE = zeros(length(subsystem_labels(trajectory)))
        TMI = 0.0
        for meas in measurements
            EE .+= meas.entropy
            TMI += meas.tmi
        end
        EE ./= trajectory.number_of_measurements
        TMI /= trajectory.number_of_measurements
        jldopen(filename, "w") do file
            file["average"] = Measurement(EE, TMI)
        end
        return Measurement(EE, TMI)
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