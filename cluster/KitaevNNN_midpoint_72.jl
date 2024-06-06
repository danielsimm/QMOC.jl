using QMOC

L = 72

pj = 1/4
pjx = pjy = pjz = pj/3
px = py = pz = 1/4
params = [px, py, pz, pjx, pjy, pjz]


n_trajectories = 768
thermalization_time = 2*L + 20
n_samples = 30
sample_distance = 2


circuits = [QMOC.KitaevNNNCircuit(L, params)]
QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevNNN_midpoint_$(L)")