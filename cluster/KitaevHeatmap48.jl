using QMOC
using JLD2
using DelimitedFiles

L = 48
params = parameter_wedge(15)
circuits = [QMOC.KitaevCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 96
thermalization_time = 3*L
n_samples = 100
sample_distance = 3


QMOC.mpi_sample_I3(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevHeatmap$(L)")
