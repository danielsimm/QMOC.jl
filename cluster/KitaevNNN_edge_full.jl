using QMOC

L = 36

pjs = unique([range(0, 0.25, length=5)..., range(0.25, 1/3, length=10)..., range(1/3, 0.5, length=10)...])
pys = pzs = zeros(length(pjs))
pxs = 1 .- pjs
pjxs = pjys = pjzs = pjs./3
params = [[pxs[i], pys[i], pzs[i], pjxs[i], pjys[i], pjzs[i]] for i in eachindex(pxs)]
circuits = [QMOC.KitaevNNNCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 192
thermalization_time = 3*L
n_samples = 100
sample_distance = 3


QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevNNN_edge_full")
