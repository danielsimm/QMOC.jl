using QMOC
using CairoMakie
using Statistics
L = 24
c = QMOC.KitaevCircuit(L, [1//3, 1//3, 1//3])

thermal_states = [QMOC.thermal_state(c) for _ in 1:10]

ent = mean([QMOC.entropy(s, c) for s in thermal_states])

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, QMOC.subsystem_labels(c), ent ./L )
fig