using QMOC
using JLD2
using DelimitedFiles

function parameter_line(pointA, pointB, steps)
    parameters = []
    step = -(pointA .- pointB)//(steps-1)
    for i in 1:steps
        px = pointA[1] + step[1]*(i-1)
        py = pointA[2] + step[2]*(i-1)
        pz = pointA[3] + step[3]*(i-1)
        pj = pointA[4] + step[4]*(i-1)
        push!(parameters, [px, py, pz, pj])
    end
    return parameters
end

L = 36
Kc = 0.672
params = parameter_line([1,0,0,0], [1//4,1//4,1//4, 1//4], 25)

push!(params, [Kc, (1-Kc)/3, (1-Kc)/3, (1-Kc)/3])
params = [[params[i][1], params[i][2], params[i][3], params[i][4]/3, params[i][4]/3, params[i][4]/3] for i in eachindex(params)]
circuits = [QMOC.KitaevNNNCircuit(L, params[i]) for i in eachindex(params)]
n_trajectories = 192
thermalization_time = 3*L
n_samples = 150
sample_distance = 3


QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "FullKitaevNNN2$(L)")
