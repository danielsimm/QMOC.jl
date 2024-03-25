module QMOC

using Dates
using LinearAlgebra
using QuantumClifford
using Statistics
using JLD2
using DelimitedFiles
using PrettyTables
using MPI
using Random

# import Distributed: @distributed, @everywhere, fetch, myid, pmap, remotecall
import JLD2: jldopen, load

include("circuits.jl")
include("Trajectories.jl")
include("Geometry.jl")
include("KitaevKekule.jl")
include("YaoKivelson.jl")
include("Chain.jl")
include("parameter_toolbox.jl")
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

# from Simulators.jl
export Simulation, simulation, simulate

# from Trajectories.jl
export run, Trajectory, Measurement

# from Analysis.jl
export evaluate

# from metadataHandling.jl
export backup, cleanup, missingTrajectories, orphanedTrajectories, archive, loadMetadata, saveMetadata, metadata

# from Parameters.jl
export parameter_full, parameter_wedge, parameter_line

end
