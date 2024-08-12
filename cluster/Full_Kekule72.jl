using QMOC
using JLD2
using DelimitedFiles


L = 72
Kc = 0.623
params = QMOC.parameter_line([0,1,0], [1//3,1//3,1//3], 25)
push!(params, [(1-Kc)/2, Kc, (1-Kc)/2])
circuits = [QMOC.KekuleCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 192
thermalization_time = 3*L
n_samples = 150
sample_distance = 3


QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "FullKitaev$(L)")
