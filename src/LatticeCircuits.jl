module LatticeCircuits

using QuantumClifford
using JLD2

include("Trajectories.jl")
include("Honeycomb.jl")
include("Decorated Honeycomb.jl")
include("Chain.jl")
include("FastProjections.jl")
include("Parameters.jl")
include("Simulators.jl")

export HClattice
# Write your package code here.

end
