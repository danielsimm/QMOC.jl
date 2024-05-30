using QMOC

L = 12

pj = 1/3
pjx = pjy = pjz = pj/3
px = 1 - pj
py = pz = 0.0
params = [px, py, pz, pjx, pjy, pjz]


n_trajectories = 384
thermalization_time = 3*L
n_samples = 100
sample_distance = 3


circuits = [QMOC.KitaevNNNCircuit(L, params)]
QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevNNN_edge_critical_$(L)")