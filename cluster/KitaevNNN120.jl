using QMOC
using JLD2
using DelimitedFiles


L = 120
Kc = 0.672
dK1 = 8.0
dK2 = 2.0
nu = 0.93
n_vals = 100 #<-- number of points in the phase diagram
K_min1 = max(0.33, Kc - dK1 / L^(1/nu))
K_min2 = max(0.33, Kc - dK2 / L^(1/nu))
K_max1 = min(1.0, Kc + dK1 / L^(1/nu))
K_max2 = min(1.0, Kc + dK2 / L^(1/nu))
rangeA = range(K_min1, K_min2, length=floor(Int, n_vals/3))
rangeB = range(K_min2, K_max2, length=ceil(Int, n_vals/3))
rangeC = range(K_max2, K_max1, length=floor(Int, n_vals/3))
pxs = reduce(vcat, [rangeA, rangeB, rangeC])
pys = pzs = pjs = (1 .- pxs)./3
pjxs = pjys = pjzs = pjs./3
params = [[pxs[i], pys[i], pzs[i], pjxs[i], pjys[i], pjzs[i]] for i in eachindex(pxs)]
circuits = [QMOC.KitaevNNNCircuit(L, params[i]) for i in eachindex(params)]

n_trajectories = 96
thermalization_time = 2*L + 20
n_samples = 1
sample_distance = 1

println("KitaevNNN L = $(L)")

QMOC.mpi_sample_I3(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, "KitaevNNN$(L)")
