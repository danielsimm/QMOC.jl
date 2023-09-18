using Distributed
addprocs(4, exeflags="-t2")
@everywhere using LatticeCircuits

simulate(LatticeCircuits.test_simulation(:YaoKivelsonXYZ))