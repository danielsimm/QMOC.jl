import Base: hash
import Graphs
import SimpleWeightedGraphs as swg

abstract type AbstractCircuit end
abstract type DecoratedHoneycombCircuit <: AbstractCircuit end
abstract type HoneycombCircuit <: AbstractCircuit end
abstract type ChainCircuit <: AbstractCircuit end

"""
    hash(c::AbstractCircuit) -> UInt

    Returns a hash of the trajectory properties.
"""
function hash(c::AbstractCircuit)
    return hash("$(typeof(c))_$(c.size)_$(c.params)")
end
struct GenericCircuit <: AbstractCircuit
    size::Int
    dims::Int
    nqubits::Int
    graph::swg.SimpleWeightedGraph
    name::String
    function GenericCircuit(size, dims, nqubits, graph, name)
        return new(size, dims, nqubits, graph, name)
    end
end

struct YaoKivelsonXYZCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonXYZCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct YaoKivelsonOrientableCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonOrientableCircuit(size, J::Real)
        return new(size, 6*size^2, [J, 1-J])
    end
    function YaoKivelsonOrientableCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct YaoKivelsonNonorientableCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonNonorientableCircuit(size, J::Real)
        return new(size, 6*size^2, [J, 1-J])
    end
    function YaoKivelsonNonorientableCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct KitaevCircuit <: HoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function KitaevCircuit(size, params)
        return new(size, 2*size^2, params)
    end
end

struct KekuleCircuit <: HoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function KekuleCircuit(size, params)
        return new(size, 2*size^2, params)
    end
end