using LatticeCircuits
using BenchmarkTools

### parameters
test_parameters = 
[
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1//3, 1//3, 1//3]
]
### trajectory
traj = LatticeCircuits.trajectory(
    :ChainPP, 
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

init = LatticeCircuits.initialise(traj)

operators = LatticeCircuits.get_operators(traj)

@benchmark LatticeCircuits.circuit!($(init), $(traj), $(operators))
