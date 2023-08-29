### Geometry ###

struct DHCSite
    lindex::Int
    xneighbour::Int
    yneighbour::Int
    zneighbour::Int
    cartesianx::Float64
    cartesiany::Float64
end

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

# function DHC_subsystem(L, ls)
#     # cut along y bonds
#     nqubits = 6*L^2
#     subsystem = []
#     for i in ls
#         A = (6*(i-1) + 1):6*L:nqubits
#         B = (6*(i-1) + 2):6*L:nqubits
#         C = (6*(i-1) + 3):6*L:nqubits
#         D = (6*(i-1) + 4):6*L:nqubits
#         E = (6*(i-1) + 5):6*L:nqubits
#         F = (6*(i-1) + 6):6*L:nqubits
#         push!(subsystem, A)
#         push!(subsystem, B)
#         push!(subsystem, C)
#         push!(subsystem, D)
#         push!(subsystem, E)
#         push!(subsystem, F)
#     end
#     return sort(union(subsystem...))
# end

function DHC_subsystem(L, ls)
    # cut along Z bonds
    sites = []
    for l in ls
        for i in ((l-1)*6*L +1):(l*6*L)
            push!(sites, i)
        end
    end
    return sites
end

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

function DHClattice(L) :: Vector{DHCSite}
    lattice = Vector{DHCSite}(undef, 6*L*L)
    for i in 1:6*L*L
        cartesian = _DHC_cartesian_position(i, L)
        lattice[i] = DHCSite(
            i,
            _DHC_xneighbour(i, L),
            _DHC_yneighbour(i, L),
            _DHC_zneighbour(i, L),
            -cartesian[1],
            cartesian[2])
    end
    return lattice
end


### Operators ###

function _DHC_largeloop_sites(L)
    loops = []
    for i in 1:6*L^2
        if mod1(i, 6) == 3
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            site4 = _DHC_yneighbour(site3, L)
            site5 = _DHC_zneighbour(site4, L)
            site6 = _DHC_xneighbour(site5, L)
            site7 = _DHC_yneighbour(site6, L)
            site8 = _DHC_zneighbour(site7, L)
            site9 = _DHC_xneighbour(site8, L)
            site10 = _DHC_yneighbour(site9, L)
            site11 = _DHC_zneighbour(site10, L)
            site12 = _DHC_xneighbour(site11, L)
            push!(loops, [site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12])
        end
    end
    return loops
end

function _DHC_smallloop_sites(L)
    loops = []
    for i in 1:6*L^2
        if mod1(i, 6) == 2
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            push!(loops, [site1, site2, site3])
        elseif mod1(i, 6) == 4
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            push!(loops, [site1, site2, site3])
        end
    end
    return loops
end

function _DHC_wilsonline_sites(L)
    loop1 = zeros(Int, 6*L)
    loop1[1] = 1
    for i in 2:6*L
        if iseven(i)
            loop1[i] = _DHC_xneighbour(loop1[i-1], L)
        else
            loop1[i] = _DHC_zneighbour(loop1[i-1], L)
        end
    end

    loop2 = zeros(Int, 6*L)
    loop2[1] = 1
    for i in 2:6*L
        if iseven(i)
            loop2[i] = _DHC_xneighbour(loop2[i-1], L)
        else
            loop2[i] = _DHC_yneighbour(loop2[i-1], L)
        end
    end

    return loop1, loop2
end

function _DHC_largeloop_operators(L) :: Vector{PauliOperator}
    loops = _DHC_largeloop_sites(L)
    N = 6*L^2
    operators = []
    for loop in loops
        Xarr = falses(N)
        Zarr = falses(N)

        Xarr[loop[1]] = true

        Xarr[loop[2]] = true
        Zarr[loop[2]] = true

        Zarr[loop[3]] = true

        Xarr[loop[4]] = true

        Xarr[loop[5]] = true
        Zarr[loop[5]] = true

        Zarr[loop[6]] = true

        Xarr[loop[7]] = true

        Xarr[loop[8]] = true
        Zarr[loop[8]] = true

        Zarr[loop[9]] = true

        Xarr[loop[10]] = true

        Xarr[loop[11]] = true
        Zarr[loop[11]] = true

        Zarr[loop[12]] = true
       push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_smallloop_operators(L) :: Vector{PauliOperator}
    loops = _DHC_smallloop_sites(L)
    N = 6*L^2
    operators = []
    for loop in loops
        Xarr = falses(N)
        Zarr = falses(N)

        Xarr[loop[1]] = true

        Xarr[loop[2]] = true
        Zarr[loop[2]] = true

        Zarr[loop[3]] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_wilsonline_operators(L) :: Vector{PauliOperator}
    loop1, loop2 = _DHC_wilsonline_sites(L)
    N = 6*L^2
    operators = Vector{PauliOperator}(undef, 2)
    Xarr = falses(N)
    Zarr = falses(N)
    # types = _generate_pattern_2(length(loop1))
    # for i in eachindex(loop1)
    #     site = loop1[i]
    #     if types[i] == :X
    #         Xarr[site] = true
    #     elseif types[i] == :Z
    #         Zarr[site] = true
    #     end
    # end
    for i in eachindex(loop1)
        site = loop1[i]
        Xarr[site] = true
        Zarr[site] = true
    end
    operators[1] = PauliOperator(0x00, Xarr, Zarr)
    Xarr = falses(N)
    Zarr = falses(N)
    # types = _generate_pattern_1(length(loop2))
    # for i in eachindex(loop2)
    #     site = loop2[i]
    #     if types[i] == :X
    #         Xarr[site] = true
    #     elseif types[i] == :Y
    #         Xarr[site] = true
    #         Zarr[site] = true
    #     end
    # end
    for i in eachindex(loop2)
        site = loop2[i]
        Zarr[site] = true
    end
    operators[2] = PauliOperator(0x00, Xarr, Zarr)
    return operators
end

function _generate_pattern_1(n)
    out = []
    patterns = [:Y, :X, :X, :Y]
    pattern_idx = 1
    for i in 1:n
        push!(out, patterns[pattern_idx])
        pattern_idx = mod(pattern_idx % 4 + 1, 5)
    end
    return out
end

function _generate_pattern_2(n)
    out = []
    patterns = [:Z, :X, :X, :Z]
    pattern_idx = 1
    for i in 1:n
        push!(out, patterns[pattern_idx])
        pattern_idx = mod(pattern_idx % 4 + 1, 5)
    end
    return out
end

_generate_pattern_1(10)

function _DHC_ZZ_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = Vector{PauliOperator}(undef, Int(N/2))
    for i in 1:2:N
        Xarr = falses(N)
        Zarr = falses(N)
        Zarr[i] = true
        Zarr[_DHC_zneighbour(i,L)] = true
        operators[Int((i+1)/2)] = PauliOperator(0x00, Xarr, Zarr)
    end
    return operators
end


# function DHC_initial_state_old(L)
#     stab = one(MixedDestabilizer, 6*L^2)
#     for op in _DHC_ZZ_operators(L)
#         project!(stab, op, keep_result=false, phases=false)
#     end
#     for op in _DHC_largeloop_operators(L)
#         project!(stab, op, keep_result=false, phases=false)
#     end
#     for op in _DHC_smallloop_operators(L)
#         project!(stab, op, keep_result=false, phases=false)
#     end
#     for op in _DHC_wilsonline_operators(L)
#         project!(stab, op, keep_result=false, phases=false)
#     end
#     if QuantumClifford.trusted_rank(stab) != 6*L^2
#         println("Error: initial state is not pure!")
#     end 
#     return stab
# end

function DHC_initial_state(L)
    stabs = [_DHC_ZZ_operators(L)..., _DHC_largeloop_operators(L)..., _DHC_smallloop_operators(L)..., _DHC_wilsonline_operators(L)...]
    return MixedDestabilizer(Stabilizer(stabs))
end

QuantumClifford.trusted_rank(DHC_initial_state(4))
nqubits(DHC_initial_state(4))

function _DHC_randomXXmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/6))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Xarr[random_site] = true
    Xarr[_DHC_xneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _DHC_randomYYmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/6))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Xarr[random_site] = true
    Zarr[random_site] = true
    Xarr[_DHC_yneighbour(random_site, L)] = true
    Zarr[_DHC_yneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _DHC_randomZZmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/6))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Zarr[random_site] = true
    Zarr[_DHC_zneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _DHC_randomsmallloopmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/6))
    loop = rand(_DHC_smallloop_sites(L))
    direction = rand([:X, :Y, :Z])
    Xarr = falses(N)
    Zarr = falses(N)

    if direction == :X
        Xarr[loop[2]] = true
        Xarr[loop[3]] = true
    elseif direction == :Y
        Xarr[loop[1]] = true
        Xarr[loop[3]] = true
        Zarr[loop[1]] = true
        Zarr[loop[3]] = true
    elseif direction == :Z
        Zarr[loop[1]] = true
        Zarr[loop[2]] = true
    end
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
end

function _DHC_randomlargeloopmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/6))
    loop = rand(_DHC_largeloop_sites(L))
    direction = rand([:X1, :X2, :Y1, :Y2, :Z1, :Z2])
    Xarr = falses(N)
    Zarr = falses(N)
    if direction == :Z1
        Zarr[loop[1]] = true
        Zarr[loop[2]] = true
    elseif direction == :Y1
        Xarr[loop[3]] = true
        Zarr[loop[3]] = true
        Xarr[loop[4]] = true
        Zarr[loop[4]] = true
    elseif direction == :X1
        Xarr[loop[5]] = true
        Xarr[loop[6]] = true
    elseif direction == :Z2
        Zarr[loop[7]] = true
        Zarr[loop[8]] = true
    elseif direction == :Y2
        Xarr[loop[9]] = true
        Zarr[loop[9]] = true
        Xarr[loop[10]] = true
        Zarr[loop[10]] = true
    elseif direction == :X2
        Xarr[loop[11]] = true
        Xarr[loop[12]] = true
    end
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
end

QuantumClifford.trusted_rank(DHC_initial_state(12))
6*12^2

### Dynamics ###

function _DHC_XYZ_timestep!(state::QuantumClifford.AbstractStabilizer, px, py, pz)
    for subtime in 1:nqubits(state)
        p = rand()
        if p < px
            _DHC_randomXXmeasurement!(state)
        elseif p < px + py
            _DHC_randomYYmeasurement!(state)
        elseif p < px + py + pz
            _DHC_randomZZmeasurement!(state)
        end
    end
    return nothing
end

function _DHC_JJ_timestep!(state::QuantumClifford.AbstractStabilizer, pJ)
    for subtime in 1:nqubits(state)
        p = rand()
        if p < pJ
            _DHC_randomsmallloopmeasurement!(state)
        else
            _DHC_randomlargeloopmeasurement!(state)
        end
    end
    return nothing
end


### Observables ###

function DHC_EE(state, L; algo=Val(:rref))
    EE = zeros(L+1)
    for i in 1:L
        EE[i+1] = entanglement_entropy(state, DHC_subsystem(L, 1:i), algo)
    end
    return EE
end

function DHC_TMI(state, L; algo=Val(:rref))
    if mod(L, 4) != 0
        @info "L must be a multiple of 4, but is $(numberOfQubits). Tripartite mutual information is ill-defined."
    end
    A = DHC_subsystem(L, 1:Int(L/4))
    B = DHC_subsystem(L, Int(L/4)+1:Int(L/2))
    C = DHC_subsystem(L, Int(L/2)+1:Int(3L/4))
    SA = entanglement_entropy(state, A, algo)
    SB = entanglement_entropy(state, B, algo)
    SC = entanglement_entropy(state, C, algo)
    SAB = entanglement_entropy(state, union(A,B), algo)
    SBC = entanglement_entropy(state, union(B,C), algo)
    SAC = entanglement_entropy(state, union(A,C), algo)
    SABC = entanglement_entropy(state, union(A, B, C), algo)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end

function DHC_observables(state; algo = Val(:rref)) # returns EE, TMI, multithreaded
    N = nqubits(state)
    L = Int(sqrt(N/6))
    if mod(L, 4) != 0
        @info "L must be a multiple of 4, but is $(N). Tripartite mutual information is ill-defined."
    end
    results = zeros(L+1+7)
    subsystems = [DHC_subsystem(L, 1:i) for i in 1:L]
    A = DHC_subsystem(L, 1:Int(L/4))
    B = DHC_subsystem(L, Int(L/4)+1:Int(L/2))
    C = DHC_subsystem(L, Int(L/2)+1:Int(3L/4))
    push!(subsystems, A)
    push!(subsystems, B)
    push!(subsystems, C)
    push!(subsystems, union(A,B))
    push!(subsystems, union(B,C))
    push!(subsystems, union(A,C))
    push!(subsystems, union(A, B, C))
    Threads.@threads for i in eachindex(subsystems)
        results[i+1] = entanglement_entropy(copy(state), subsystems[i], algo)
    end
    EE = results[1:L+1]
    TMI = results[L+2] + results[L+3] + results[L+4] - results[L+5] - results[L+6] - results[L+7] + results[L+8]
    return EE, TMI
end


### Benchmark
state = DHC_initial_state(48)
function DHC_benchmark_evolution(state)
    for i in 1:3
        _DHC_XYZ_timestep!(state, 0.1, 0.1, 0.1)
    end
    return nothing
end

function DHC_benchmark_observables(state)
    EE, TMI = DHC_observables(state)
    return nothing
end

function DHC_benchmark_async(state)
    for _ in 1:10
        for i in 1:3
            _DHC_XYZ_timestep!(state, 0.1, 0.1, 0.1)
        end
        @async DHC_observables(state)
    end
end

function DHC_benchmark_serial(state)
    for _ in 1:10
        for i in 1:3
            _DHC_XYZ_timestep!(state, 0.1, 0.1, 0.1)
        end
        DHC_observables(state)
    end
end

@benchmark DHC_benchmark_async(DHC_initial_state(12))

@benchmark DHC_benchmark_serial(DHC_initial_state(12))