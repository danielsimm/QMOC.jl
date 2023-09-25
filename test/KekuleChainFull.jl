using Distributed
addprocs(6, exeflags="-t2")
@everywhere using LatticeCircuits
parameters = parameter_full(10)
simulate(simulation(:ChainKekule, 120, "KekuleChainFull128", 240, parameters, checkpoints=false, verbosity=:none))
