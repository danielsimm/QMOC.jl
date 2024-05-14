using QMOC

L = 32
Kc = 0.66
dK = 8.0
nu = 1.0
n_vals = 100 #<-- number of points in the phase diagram
K_min = max(0.33, Kc - dK / L^(1/nu))
K_max = min(1.0, Kc + dK / L^(1/nu))
pxs = range(K_min, K_max, length=n_vals)
pys = pzs = (1 .- pxs)./2
params = [[pxs[i], pys[i], pzs[i]] for i in eachindex(pxs)]
circuits = [QMOC.KitaevCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 192
thermalization_time = 3*L
n_samples = 150
sample_distance = 3


QMOC.mpi_sample_I3(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "Kitaev$(L)")
