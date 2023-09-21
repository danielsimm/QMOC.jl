using Distributed
addprocs(6, exeflags="-t2")
@everywhere using LatticeCircuits
parameters = parameter_line([1//2, 1//2, 0], :center, 40)
PPline256 = simulation(:ChainPP, 256, "PPcritical256", 240, parameters, checkpoints=true)
PPline512 = simulation(:ChainPP, 512, "PPcritical512", 240, parameters, checkpoints=true)
PPline1024 = simulation(:ChainPP, 1024, "PPcritical1024", 240, parameters, checkpoints=true)

simulate(PPline256)
simulate(PPline512)
simulate(PPline1024)