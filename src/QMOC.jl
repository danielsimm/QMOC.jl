module QMOC

using Dates
using LinearAlgebra
using QuantumClifford
using Statistics
using JLD2
using PrettyTables

import Distributed: @distributed, @everywhere, fetch, myid, pmap, remotecall
import JLD2: jldopen, load


include("Trajectories.jl")
include("Geometry.jl")
include("KitaevKekule.jl")
include("YaoKivelson.jl")
include("Chain.jl")
include("Parameters.jl")
include("Simulators.jl")
include("Analysis.jl")

if !isdir("data")
    mkdir("data")
end
if !isdir("data/metadata")
    mkdir("data/metadata")
end
include("MetadataHandling.jl")
# checkMetadataIntegrity()

println("QMOC.jl loaded on worker $(myid()) with $(Threads.nthreads()) threads.")

# from Simulators.jl
export Simulation, simulation, simulate

# from Trajectories.jl
export run, Trajectory, Measurement

# from Analysis.jl
export evaluate

# from metadataHandling.jl
export metadata, loadSimulation, missingTrajectories

# from Parameters.jl
export parameter_full, parameter_wedge, parameter_line

metadata()

end
