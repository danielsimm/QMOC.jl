using Distributed
addprocs()
@everywhere include("src/LatticeCircuits.jl")
@everywhere using .LatticeCircuits
sim = simulation(:Kekule, 36, "test", 16, [[1, 0, 0]], checkpoints=true, verbosity=:high)
simulate(sim)