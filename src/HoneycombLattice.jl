function _HC_xneighbour(i::Int64, L::Int64)
    row = div(i-1, 2*L)
    if iseven(i)
        return mod1(i-1, 2*L) + row*2*L
    else
        return mod1(i+1, 2*L) + row*2*L
    end
end

function _HC_yneighbour(i::Int64, L::Int64)
    row = div(i-1, 2*L)
    if iseven(i)
        return mod1(i+1, 2*L) + row*2*L
    else
        return mod1(i-1, 2*L) + row*2*L
    end
end

function _HC_zneighbour(i::Int64, L::Int64)
    if iseven(i)
        return mod1(i+2*L-1, 2*L^2)
    else
        return mod1(i-2*L+1, 2*L^2)
    end
end

function _HC_neighbour(i::Int64, L::Int64, dir::Symbol)
    if dir == :X
        return _HC_xneighbour(i, L)
    elseif dir == :Y
        return _HC_yneighbour(i, L)
    elseif dir == :Z
        return _HC_zneighbour(i, L)
    end
end

function _HC_redneighbour(i::Int64, L::Int64)
    dir = _HC_kekuledirection(i, L, :red)
    return _HC_neighbour(i, L, dir), dir
end

function _HC_greenneighbour(i::Int64, L::Int64)
    dir = _HC_kekuledirection(i, L, :green)
    return _HC_neighbour(i, L, dir), dir
end

function _HC_blueneighbour(i::Int64, L::Int64)
    dir =  _HC_kekuledirection(i, L, :blue)
    return _HC_neighbour(i, L, dir), dir
end

# function HC_subsystem(L,ls)
#     qubits = []
#     for l in ls
#         push!(qubits, collect((2L*(l-1)+1):2L*l))
#     end
#     return reduce(vcat, qubits)
# end

function _HC_xyindex(i::Int64, L::Int64)
    y = div(i-1, 2*L) + 1
    x = mod1(i, 2*L)
    return x, y
end

function _HC_kekuledirection(i::Int64, L::Int64, color::Symbol)
    redevenmatrix = [:Z :X :Y 
                  :Y :Z :X 
                  :X :Y :Z]
    redoddmatrix = [:Y :X :Z 
                 :Z :Y :X 
                 :X :Z :Y]
    greenevenmatrix = [:Y :Z :X 
                    :X :Y :Z 
                    :Z :X :Y]
    greenoddmatrix = [:Z :Y :X
                     :X :Z :Y
                     :Y :X :Z]
    blueevenmatrix = [:X :Y :Z
                    :Z :X :Y
                    :Y :Z :X]
    blueoddmatrix = [:X :Z :Y
                   :Y :X :Z
                   :Z :Y :X]
    x, y = _HC_xyindex(i, L)
    y = mod1(y, 3)
    
    if color == :red
        if iseven(x)
            x = Int64(mod1(x, 6)/2)
            return redevenmatrix[y,x]
        else
            x = Int64(mod1(x+1, 6)/2)
            return redoddmatrix[y,x]
        end
    elseif color == :green
        if iseven(x)
            x = Int64(mod1(x, 6)/2)
            return greenevenmatrix[y,x]
        else
            x = Int64(mod1(x+1, 6)/2)
            return greenoddmatrix[y,x]
        end
    elseif color == :blue
        if iseven(x)
            x = Int64(mod1(x, 6)/2)
            return blueevenmatrix[y,x]
        else
            x = Int64(mod1(x+1, 6)/2)
            return blueoddmatrix[y,x]
        end
    end
end

function HC_subsystem(L, ls)
    qubits = []
    for l in ls
        start = 2*(l-1)+1
        for i in 1:2*L
            push!(qubits, start)
            if iseven(start)
                start = _HC_zneighbour(start, L)
            else
                start = _HC_xneighbour(start, L)
            end
        end
    end
    return sort(qubits)
end

function _HC_cartesian_position(i::Int64, L::Int64)
    if i == 1
        return [0, 0]
    end
    row = div(i-1, 2*L)
    if iseven(i)
        return _HC_cartesian_position(i-1, L) + [1, sqrt(3)/2]
    elseif row > 0
        return _HC_cartesian_position(_HC_zneighbour(i, L), L) + [0, sqrt(7/4)]
    else
        return _HC_cartesian_position(_HC_yneighbour(i, L), L) + [1, -sqrt(3)/2]
    end
end


struct HoneycombLatticeSite
    index::Int64
    xneighbour::Int64
    yneighbour::Int64
    zneighbour::Int64
    redneighbour::Int64
    greenneighbour::Int64
    blueneighbour::Int64
    cartesianx::Float64
    cartesiany::Float64
end

struct HoneycombLattice <: Lattice
    L::Int
    nqubits::Int
    sites::Vector{HoneycombLatticeSite}
end

function honeycomblattice(L) :: HoneycombLattice
    sites = Vector{HoneycombLatticeSite}(undef, 2*L^2)
    for i in eachindex(sites)
        sites[i] = HoneycombLatticeSite(
            i, 
            _HC_xneighbour(i, L), 
            _HC_yneighbour(i, L), 
            _HC_zneighbour(i, L), 
            _HC_redneighbour(i, L)[1], 
            _HC_greenneighbour(i, L)[1], 
            _HC_blueneighbour(i, L)[1], 
            _HC_cartesian_position(i, L)[1],
            _HC_cartesian_position(i, L)[2]
        )
    end
    return HoneycombLattice(L, 2*L^2, sites)
end