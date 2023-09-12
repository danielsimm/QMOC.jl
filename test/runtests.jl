using LatticeCircuits

### parameters
test_parameters = 
[
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1//3, 1//3, 1//3]
]
### single trajectory test
traj = LatticeCircuits.trajectory(
    :ChainPPFast, 
    1024,
    1024, 
    "trajectory_test", 
    test_parameters[4], 
    true, 
    :debug, 
    1, 
    3*1024, 
    30, 
    10)

run(traj)

# using Distributed
# addprocs()
# # @everywhere include("src/LatticeCircuits.jl")
# @everywhere using LatticeCircuits
# sim = simulation(:Kekule, 36, "test", 16, [[1, 0, 0]], checkpoints=true, verbosity=:high)
# simulate(sim)