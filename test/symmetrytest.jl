include("LatticeCircuits.jl")
using CairoMakie
include("src/Plotting.jl")


function thermalize!(state, L, px, py, pz)
    for time in 1:6*L
        _DHC_XYZ_timestep!(state, px, py, pz)
    end
end

function _yaokivelson_XYZ_trajectory_test(init, L, parameters)
    px = parameters[1]
    py = parameters[2]
    pz = parameters[3]
    measurements = 100

    state = init
    for layer in 1:6*L
        _DHC_XYZ_timestep!(state, px, py, pz)
    end

    ### Measurements
    entropy = zeros(measurements, L+1)
    TMI = zeros(measurements)
    for measurement in 1:measurements
        for layer in 1:3
            _DHC_XYZ_timestep!(state, px, py, pz)
        end
        entropy[measurement, :] = DHC_EE(state, L)
        TMI[measurement] = DHC_TMI(state, L)
    end

    return sum(entropy, dims=1)./measurements, sum(TMI)./measurements
end

function generate_phase_diagram_test(L::Int64, initial_direction::Symbol)
    # Generate a phase diagram for a given L and initial direction
    # initial_direction can be :X, :Y, or :Z
    # Returns a figure
    params = parameter_full(15)
    TMIs = zeros(length(params))
    Threads.@threads for i in eachindex(params)
        TMIs[i] = _yaokivelson_XYZ_trajectory_test(DHC_initial_state(L, initial_direction), L, params[i])[2]
    end
    return voronoi_tesselation_plot(parametric_to_cartesian.(params), TMIs)
end

generate_phase_diagram_test(16, :Z)



