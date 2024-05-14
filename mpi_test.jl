using QMOC

L = 12
Kc = 0.66
dK = 3.0
nu = 1.0
n_vals = 100 #<-- number of points in the phase diagram
K_min = max(0, Kc - dK / L^(1/nu))
K_max = Kc + 0.3* dK / L^(1/nu)
pxs = range(K_min, K_max, length=n_vals)
pys = pzs = (1 .- pxs)./2
params = [[pxs[i], pys[i], pzs[i]] for i in eachindex(pxs)]
circuits = [QMOC.KitaevCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 128
thermalization_time = 3*L
n_samples = 100
sample_distance = 3


QMOC.mpi_sample_I3(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "test")
