abstract type HoneycombTrajectory <: Trajectory end
struct KitaevTrajectory <: HoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Any}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

struct KekuleTrajectory <: HoneycombTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Any}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end


### Operators ###

function _HC_ZZ_operator(i::Int64, L::Int64)
    j = _HC_zneighbour(i, L)
    Xarr = falses(2*L^2)
    Zarr = falses(2*L^2)
    Zarr[i] = true
    Zarr[j] = true
    return PauliOperator(0x00, Xarr, Zarr)
end

function _HC_ZZ_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2-1)
    for i in eachindex(operators)
        operators[i] = _HC_ZZ_operator(2*i, L)
    end
    return operators
end

function _HC_WilsonPlaquette_operator(i::Int64, L::Int64)
    i = i #X
    j = _HC_zneighbour(i, L) #Y
    k = _HC_xneighbour(j, L) #Z
    l = _HC_yneighbour(k, L) #X
    m = _HC_zneighbour(l, L) #Y
    n = _HC_xneighbour(m, L) #Z
    Xarr = falses(2*L^2)
    Zarr = falses(2*L^2)

    Xarr[i] = true

    Xarr[j] = true
    Zarr[j] = true

    Zarr[k] = true

    Xarr[l] = true

    Xarr[m] = true
    Zarr[m] = true

    Zarr[n] = true

    return PauliOperator(0x00, Xarr, Zarr)
end

function _HC_WilsonPlaquette_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2-1)
    for i in eachindex(operators)
        operators[i] = _HC_WilsonPlaquette_operator(2*i, L)
    end
    return operators
end

# function XZWilsonLines(L::Int64)
#     Lines = []
#     for i in 1:L
#         site = 2*i-1
#         Xarr = falses(2*L^2)
#         Zarr = falses(2*L^2)
#         for j in 1:2*L
#             Xarr[site] = true
#             Zarr[site] = true
#             if iseven(j)
#                 site = _HC_zneighbour(site, L)
#             else
#                 site = _HC_xneighbour(site, L)
#             end
#         end
#         push!(Lines, PauliOperator(0x00, Xarr, Zarr))
#     end
#     return Lines
# end

# function XYWilsonLines(L::Int64)
#     Lines = []
#     start = 1
#     for i in 1:L
#         site = start
#         Xarr = falses(2*L^2)
#         Zarr = falses(2*L^2)
#         for j in 1:2*L
#             Zarr[site] = true
#             if iseven(j)
#                 site = _HC_yneighbour(site, L)
#             else
#                 site = _HC_xneighbour(site, L)
#             end
#         end
#         push!(Lines, PauliOperator(0x00, Xarr, Zarr))
#         start = _HC_xneighbour(start, L)
#         start = _HC_zneighbour(start, L)
#     end
#     return Lines
# end

function _HC_WilsonLoops(L::Int64)
    site1 = 1
    site2 = 1
    Xarr1 = falses(2*L^2)
    Zarr1 = falses(2*L^2)
    Xarr2 = falses(2*L^2)
    Zarr2 = falses(2*L^2)
    for i in 1:2*L
        Xarr1[site1] = true #XXXX...
        Zarr2[site2] = true #ZZZZ...
        if isodd(i)
            site1 = _HC_zneighbour(site1, L)
            site2 = _HC_xneighbour(site2, L)
        else
            site1 = _HC_yneighbour(site1, L)
            site2 = _HC_yneighbour(site2, L)
        end
    end
    return [PauliOperator(0x00, Xarr1, Zarr1) PauliOperator(0x00, Xarr2, Zarr2)]
end

# function HC_initial_state(L::Int64)
#     numberOfQubits = 2*L^2
#     state = one(MixedDestabilizer, numberOfQubits)
#     # for op in _HC_ZZ_operators(L)
#     #     project!(state, op, keep_result=false, phases=false)
#     # end
#     for op in _HC_WilsonPlaquette_operators(L)
#         project!(state, op, keep_result=false, phases=false)
#     end
#     for op in _HC_WilsonLoops(L)
#         project!(state, op, keep_result=false, phases=false)
#     end
#     if QuantumClifford.trusted_rank(state) != numberOfQubits
#         println("Error: initial state is not pure!")
#     end
#     return MixedDestabilizer(state)
# end

function HC_initial_state(L::Int64)
    stabs = [_HC_ZZ_operators(L)..., _HC_WilsonPlaquette_operators(L)..., _HC_WilsonLoops(L)...]
    return MixedDestabilizer(Stabilizer(stabs))
end

function initialise(trajectory::HoneycombTrajectory)
    return HC_initial_state(trajectory.size)
end

function _HC_randomXXmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Xarr[random_site] = true
    Xarr[_HC_xneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _HC_randomYYmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Zarr[random_site] = true
    Zarr[_HC_zneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _HC_randomZZmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    Xarr[random_site] = true
    Zarr[random_site] = true
    Xarr[_HC_yneighbour(random_site, L)] = true
    Zarr[_HC_yneighbour(random_site, L)] = true
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _HC_randomredmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    neighbour, direction = _HC_redneighbour(random_site, L)
    
    if direction == :X
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Y
        Zarr[random_site] = true
        Zarr[neighbour] = true
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Z
        Zarr[random_site] = true
        Zarr[neighbour] = true
    end
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _HC_randomgreenmeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    neighbour, direction = _HC_greenneighbour(random_site, L)
    
    if direction == :X
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Y
        Zarr[random_site] = true
        Zarr[neighbour] = true
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Z
        Zarr[random_site] = true
        Zarr[neighbour] = true
    end
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end

function _HC_randombluemeasurement!(state::QuantumClifford.AbstractStabilizer)
    N = nqubits(state)
    L = Int(sqrt(N/2))
    Xarr = falses(N)
    Zarr = falses(N)
    random_site = rand(1:N)
    neighbour, direction = _HC_blueneighbour(random_site, L)
    
    if direction == :X
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Y
        Zarr[random_site] = true
        Zarr[neighbour] = true
        Xarr[random_site] = true
        Xarr[neighbour] = true
    elseif direction == :Z
        Zarr[random_site] = true
        Zarr[neighbour] = true
    end
    project!(state, PauliOperator(0x00, Xarr, Zarr), keep_result=false, phases=false)
    return nothing
end


### Dynamics ###

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::KitaevTrajectory)
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < px
            _HC_randomXXmeasurement!(state)
        elseif p < px + py
            _HC_randomYYmeasurement!(state)
        elseif p < px + py + pz
            _HC_randomZZmeasurement!(state)
        end
    end
    return nothing
end

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::KekuleTrajectory)
    pr = trajectory.params[1]
    pg = trajectory.params[2]
    pb = trajectory.params[3]
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < pr
            _HC_randomredmeasurement!(state)
        elseif p < pr + pg
            _HC_randomgreenmeasurement!(state)
        else
            _HC_randombluemeasurement!(state)
        end
    end
    return nothing
end


### Observables ###

function entropy(state::QuantumClifford.MixedDestabilizer, trajectory::HoneycombTrajectory; algo=Val(:rref))
    L = trajectory.size
    EE = zeros(L+1)
    for i in 1:L
        EE[i+1] = entanglement_entropy(state, HC_subsystem(L, 1:i), algo)
    end
    return EE
end

function tmi(state::QuantumClifford.MixedDestabilizer, trajectory::HoneycombTrajectory; algo=Val(:rref))
    L = trajectory.size
    if mod(L, 4) != 0
        @info "L must be a multiple of 4, but is $(L). Tripartite mutual information is ill-defined."
    end
    A = HC_subsystem(L, 1:Int(L/4))
    B = HC_subsystem(L, Int(L/4)+1:Int(L/2))
    C = HC_subsystem(L, Int(L/2)+1:Int(3L/4))
    SA = entanglement_entropy(state, A, algo)
    SB = entanglement_entropy(state, B, algo)
    SC = entanglement_entropy(state, C, algo)
    SAB = entanglement_entropy(state, union(A,B), algo)
    SBC = entanglement_entropy(state, union(B,C), algo)
    SAC = entanglement_entropy(state, union(A,C), algo)
    SABC = entanglement_entropy(state, union(A, B, C), algo)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end