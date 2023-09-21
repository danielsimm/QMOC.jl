using Distributed
addprocs(6, exeflags="-t2")
@everywhere using LatticeCircuits
parameters = parameter_line([6//10, 4//10, 0], [4//10, 6//10, 0], 40)
PPline128 = simulation(:ChainPP, 128, "PPphasetransition256", 400, parameters, checkpoints=false)
PPline256 = simulation(:ChainPP, 256, "PPphasetransition256", 400, parameters, checkpoints=false)
PPline512 = simulation(:ChainPP, 512, "PPphasetransition512", 400, parameters, checkpoints=false)

simulate(PPline128)
simulate(PPline256)
simulate(PPline512)