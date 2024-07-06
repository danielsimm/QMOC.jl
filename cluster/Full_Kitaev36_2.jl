using QMOC
using JLD2
using DelimitedFiles


L = 36
Kc = 0.6524
params = QMOC.parameter_line([1,0,0], [1//3,1//3,1//3], 25)
push!(params, [Kc, (1-Kc)/2, (1-Kc)/2])
circuits = [QMOC.KitaevCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 192
thermalization_time = 3*L
n_samples = 150
sample_distance = 3


QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "FullKitaev2$(L)")
