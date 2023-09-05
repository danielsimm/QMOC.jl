using QuantumClifford
includet("src/Honeycomb.jl")
includet("src/Decorated Honeycomb.jl")

function _is_selfconsistent(lattice::Vector{HCSite})
    for i in eachindex(lattice)
        if i != lattice[lattice[i].xneighbour].xneighbour
            return false, i
        end
        if i != lattice[lattice[i].yneighbour].yneighbour
            return false, i
        end
        if i != lattice[lattice[i].zneighbour].zneighbour
            return false, i
        end
        if i != lattice[lattice[i].redneighbour].redneighbour
            return false, i
        end
        if i != lattice[lattice[i].blueneighbour].blueneighbour
            return false, i
        end
        if i != lattice[lattice[i].greenneighbour].greenneighbour
            return false, i
        end
    end
    return true, nothing
end

function _is_selfconsistent(lattice::Vector{DHCSite})
    for i in eachindex(lattice)
        if i != lattice[lattice[i].xneighbour].xneighbour
            return false, i
        end
        if i != lattice[lattice[i].yneighbour].yneighbour
            return false, i
        end
        if i != lattice[lattice[i].zneighbour].zneighbour
            return false, i
        end
    end
    return true, nothing
end

_is_selfconsistent(DHClattice(4))