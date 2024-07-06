using QMOC
using JLD2
using DelimitedFiles

L = 36
Kc = 0.6524
dK1 = 8.0
dK2 = 2.0
nu = 1.0
n_vals = 25 #<-- number of points in the phase diagram
K_min1 = max(0.33, Kc - dK1 / L^(1/nu))
K_min2 = max(0.33, Kc - dK2 / L^(1/nu))
K_max1 = min(1.0, Kc + dK1 / L^(1/nu))
K_max2 = min(1.0, Kc + dK2 / L^(1/nu))
rangeA = range(K_min1, K_min2, length=floor(Int, n_vals/3))
rangeB = range(K_min2, K_max2, length=ceil(Int, n_vals/3))
rangeC = range(K_max2, K_max1, length=floor(Int, n_vals/3))
pxs = reduce(vcat, [rangeA, rangeB, rangeC])
push!(pxs, 1/3)
pys = pzs = (1 .- pxs)./2
params = [[pxs[i], pys[i], pzs[i]] for i in eachindex(pxs)]
circuits = [QMOC.KitaevCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 192
thermalization_time = 3*L
n_samples = 150
sample_distance = 3


QMOC.mpi_sample_full(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "FullKitaev$(L)")
