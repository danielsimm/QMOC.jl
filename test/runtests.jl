using Distributed
addprocs(4, exeflags="-t2")
@everywhere using LatticeCircuits

simulate(LatticeCircuits.test_simulation(:Kitaev))
simulate(LatticeCircuits.test_simulation(:Kekule))
simulate(LatticeCircuits.test_simulation(:YaoKivelsonXYZ))
simulate(LatticeCircuits.test_simulation(:YaoKivelsonJJ))

printMetadata()

LatticeCircuits.loadSimulation("test_YaoKivelsonXYZ")

LatticeCircuits.test_trajectory(:Kitaev)