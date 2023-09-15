module LatticeCircuits

using Dates
using LinearAlgebra
using QuantumClifford
using Statistics
using JLD2
using PrettyTables

import Distributed: @distributed, @everywhere, fetch, myid, pmap, remotecall
import JLD2: jldopen, load

const PLOTTING = false
global METADATA_INTEGRITY = true ::Bool

include("Trajectories.jl")
include("Geometry.jl")
include("KitaevKekule.jl")
include("YaoKivelson.jl")
include("Chain.jl")
include("Parameters.jl")
include("Simulators.jl")
if PLOTTING
    include("Plotting.jl")
end

include("metadataHandling.jl")
checkMetadataIntegrity()
SimulationArchive = loadSimulationArchive()

@info "LatticeCircuits.jl loaded on worker $(myid()) with $(Threads.nthreads()) threads."

# from Simulators.jl
export Simulation, simulation, simulate

# from Trajectories.jl
export run, Trajectory, Measurement

# from metadataHandling.jl
export printMetadata

export SimulationArchive

end
