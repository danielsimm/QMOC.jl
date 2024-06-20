using QMOC

L = 40

pj = 1/2
pjx = pjy = pjz = pj/3
px = py = pz = 1/6
params = [px, py, pz, pjx, pjy, pjz]


n_trajectories = 768
thermalization_time = 3*L
n_samples = 30
sample_distance = 2


circuits = [QMOC.KitaevNNNCircuit(L, params)]
QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevNNN_highpoint_$(L)")