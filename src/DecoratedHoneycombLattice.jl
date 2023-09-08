function _DHC_xneighbour(globalindex, L)
    localindex = mod1(globalindex, 6)
    if localindex == 1
        return globalindex + 2
    elseif localindex == 2
        return globalindex + 2
    elseif localindex == 3
        return globalindex - 2
    elseif localindex == 4
        return globalindex - 2
    elseif localindex == 5
        return globalindex + 1
    elseif localindex == 6
        return globalindex - 1
    end
end

function _DHC_yneighbour(globalindex, L)
    localindex = mod1(globalindex, 6)
    if localindex == 1
        if globalindex-2 < div(globalindex, 6*L)*6*L
            return globalindex - 2 + 6*L
        else
            return globalindex - 2
        end
    elseif localindex == 2
        return globalindex + 1
    elseif localindex == 3
        return globalindex - 1
    elseif localindex == 4
        return globalindex + 2
    elseif localindex == 5
        if globalindex+2 > (div(globalindex, 6*L)+1)*6*L
            return globalindex + 2 - 6*L
        else
            return globalindex + 2
        end
    elseif localindex == 6
        return globalindex - 2
    end
end

function _DHC_zneighbour(globalindex, L)
    localindex = mod1(globalindex, 6)
    if localindex == 1
        return globalindex + 1
    elseif localindex == 2
        return globalindex - 1
    elseif localindex == 3
        return mod1(globalindex + 3 + 6*L, 6*L^2)
    elseif localindex == 4
        return globalindex + 1
    elseif localindex == 5
        return globalindex - 1
    elseif localindex == 6
        return mod1(globalindex - 3 - 6*L, 6*L^2)
    end
end

function DHC_subsystem(L, ls)
    # cut along y bonds
    nqubits = 6*L^2
    subsystem = []
    for i in ls
        A = (6*(i-1) + 1):6*L:nqubits
        B = (6*(i-1) + 2):6*L:nqubits
        C = (6*(i-1) + 3):6*L:nqubits
        D = (6*(i-1) + 4):6*L:nqubits
        E = (6*(i-1) + 5):6*L:nqubits
        F = (6*(i-1) + 6):6*L:nqubits
        push!(subsystem, A)
        push!(subsystem, B)
        push!(subsystem, C)
        push!(subsystem, D)
        push!(subsystem, E)
        push!(subsystem, F)
    end
    return sort(union(subsystem...))
end

# function DHC_subsystem(L, ls)
#     # cut along Z bonds
#     sites = []
#     for l in ls
#         for i in ((l-1)*6*L +1):(l*6*L)
#             push!(sites, i)
#         end
#     end
#     return sites
# end

function _DHC_cartesian_position(i::Int64, L::Int64)
    if i > 6
        cell = div(i-1, 6)+1
        xshift = (mod1(cell, L)-1) * [3, 0]
        row = div(i-1, 6*L)
        rowshift = row .* [-1.5, 1.5*sqrt(3)+1]
        return _DHC_cartesian_position(mod1(i, 6), L) + xshift + rowshift
    else
        if i == 1
            return [0.0, 0.0]
        elseif i == 2
            return [1.0, 0.0]
        elseif i == 3
            return [0.5, sqrt(3)/2]
        elseif i == 4
            return [1.5, -sqrt(3)/2]
        elseif i == 5
            return [2.5, -sqrt(3)/2]
        elseif i == 6
            return [2.0, -sqrt(3)]
        end
    end
end

struct DecoratedHoneycombLatticeSite
    index::Int
    xneighbour::Int
    yneighbour::Int
    zneighbour::Int
    cartesianx::Float64
    cartesiany::Float64
end

struct DecoratedHoneycombLattice <: Lattice
    L::Int
    nqubits::Int
    sites::Vector{DecoratedHoneycombLatticeSite}
end

function decoratedhoneycomblattice(L) :: DecoratedHoneycombLattice
    sites = Vector{DecoratedHoneycombLatticeSite}(undef, 6*L^2)
    for i in 1:6*L^2
        cartesian = _DHC_cartesian_position(i, L)
        sites[i] = DecoratedHoneycombLatticeSite(
            i,
            _DHC_xneighbour(i, L),
            _DHC_yneighbour(i, L),
            _DHC_zneighbour(i, L),
            -cartesian[1],
            cartesian[2])
    end
    return DecoratedHoneycombLattice(L, 6*L^2, sites)
end
