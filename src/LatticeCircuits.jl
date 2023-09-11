module LatticeCircuits

using QuantumClifford
using JLD2
using Statistics
using Dates
using LinearAlgebra
import Distributed: @distributed, pmap, remotecall, fetch, @everywhere

include("Trajectories.jl")
include("Geometry.jl")
include("KitaevKekule.jl")
include("YaoKivelson.jl")
include("Chain.jl")
include("FastProjections.jl")
include("Parameters.jl")
include("Simulators.jl")

include("metadataHandling.jl")
checkMetadataIntegrity()

# from Simulators.jl
export Simulation, simulation, simulate

# from Trajectories.jl
export run, Trajectory, Measurement

# from metadataHandling.jl
export printMetadata

end
