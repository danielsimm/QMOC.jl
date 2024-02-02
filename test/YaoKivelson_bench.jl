using QMOC
using BenchmarkTools

traj = QMOC.YaoKivelsonOrientableTest(12)

new_init = QMOC.initialise(traj)
old_init = QMOC.initialise_mixed(traj)

old_init == new_init