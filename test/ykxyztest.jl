using Distributed
addprocs(4, exeflags="-t2")
@everywhere using LatticeCircuits

sim = simulation(:YaoKivelsonXYZ, 36, "YaoKivelsonIsotropic36", 100, [[1//3, 1//3, 1//3]], checkpoints=true)

simulate(sim)