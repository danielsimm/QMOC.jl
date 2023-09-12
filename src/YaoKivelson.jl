abstract type DecoratedHoneycombTrajectory <: Trajectory end
struct YaoKivelsonXYZTrajectory <: DecoratedHoneycombTrajectory
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
    XX_projectors::Vector{PauliOperator}
    YY_projectors::Vector{PauliOperator}
    ZZ_projectors::Vector{PauliOperator}
end

function _YaoKivelsonXYZTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements)
    XX_projectors = _DHC_XX_operators(size)
    YY_projectors = _DHC_YY_operators(size)
    ZZ_projectors = _DHC_ZZ_operators(size)
    return YaoKivelsonXYZTrajectory(size, nqubits, name, parameters, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements, XX_projectors, YY_projectors, ZZ_projectors)
end

struct YaoKivelsonJJTrajectory <: DecoratedHoneycombTrajectory
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

export YaoKivelsonXYZTrajectory, YaoKivelsonJJTrajectory, DecoratedHoneycombTrajectory

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
    for i in eachindex(loop1)
        site = loop1[i]
        Xarr[site] = true
        Zarr[site] = true
    end
    operators[1] = PauliOperator(0x00, Xarr, Zarr)
    Xarr = falses(N)
    Zarr = falses(N)
    for i in eachindex(loop2)
        site = loop2[i]
        Zarr[site] = true
    end
    operators[2] = PauliOperator(0x00, Xarr, Zarr)
    return operators
end

function _DHC_ZZ_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:2:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Zarr[i] = true
        Zarr[_DHC_zneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_XX_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Xarr[i] = true
        Xarr[_DHC_xneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return unique(operators)
end

function _DHC_YY_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Xarr[i] = true
        Xarr[_DHC_yneighbour(i,L)] = true
        Zarr[i] = true
        Zarr[_DHC_yneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return unique(operators)
end

function initialise(trajectory::DecoratedHoneycombTrajectory; basis=:Z) :: MixedDestabilizer
    L = trajectory.size
    if basis == :Z
        bilinears = _DHC_ZZ_operators(L)
    elseif basis == :X
        bilinears = _DHC_XX_operators(L)
    elseif basis == :Y
        bilinears = _DHC_YY_operators(L)
    end
    stabs = [bilinears..., _DHC_largeloop_operators(L)..., _DHC_smallloop_operators(L)..., _DHC_wilsonline_operators(L)...]
    state = MixedDestabilizer(Stabilizer(stabs))
    if QuantumClifford.trusted_rank(state) != trajectory.nqubits
        @warn "Initial state is not pure."
    end
    return state
end

function _DHC_randomXXmeasurement!(state::QuantumClifford.MixedDestabilizer)
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

function _DHC_randomYYmeasurement!(state::QuantumClifford.MixedDestabilizer)
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

function _DHC_randomZZmeasurement!(state::QuantumClifford.MixedDestabilizer)
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

function _DHC_randomsmallloopmeasurement!(state::QuantumClifford.MixedDestabilizer)
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

function _DHC_randomlargeloopmeasurement!(state::QuantumClifford.MixedDestabilizer)
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


### Dynamics ###

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonXYZTrajectory)
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < px
            project!(state, rand(trajectory.XX_projectors), keep_result=false, phases=false)
        elseif p < px + py
            project!(state, rand(trajectory.YY_projectors), keep_result=false, phases=false)
        elseif p < px + py + pz
            project!(state, rand(trajectory.ZZ_projectors), keep_result=false, phases=false)
        end
    end
    return nothing
end

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonJJTrajectory)
    pJ = trajectory.params[1]
    for subtime in 1:trajectory.nqubits
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

function entropy(state::QuantumClifford.MixedDestabilizer, trajectory::DecoratedHoneycombTrajectory)
    algo=Val(:rref)
    L = trajectory.size
    EE = zeros(L+1)
    for i in 1:L
        EE[i+1] = entanglement_entropy(state, DHC_subsystem(L, 1:i), algo)
    end
    return EE
end

function tmi(state::QuantumClifford.MixedDestabilizer, trajectory::DecoratedHoneycombTrajectory)
    algo=Val(:rref)
    L = trajectory.size
    if mod(L, 4) != 0
        @info "L must be a multiple of 4, but is $(L). Tripartite mutual information is ill-defined."
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

# function DHC_observables(state; algo = Val(:rref)) # returns EE, TMI, multithreaded
#     N = nqubits(state)
#     L = Int(sqrt(N/6))
#     if mod(L, 4) != 0
#         @info "L must be a multiple of 4, but is $(N). Tripartite mutual information is ill-defined."
#     end
#     results = zeros(L+1+7)
#     subsystems = [DHC_subsystem(L, 1:i) for i in 1:L]
#     A = DHC_subsystem(L, 1:Int(L/4))
#     B = DHC_subsystem(L, Int(L/4)+1:Int(L/2))
#     C = DHC_subsystem(L, Int(L/2)+1:Int(3L/4))
#     push!(subsystems, A)
#     push!(subsystems, B)
#     push!(subsystems, C)
#     push!(subsystems, union(A,B))
#     push!(subsystems, union(B,C))
#     push!(subsystems, union(A,C))
#     push!(subsystems, union(A, B, C))
#     Threads.@threads for i in eachindex(subsystems)
#         results[i+1] = entanglement_entropy(copy(state), subsystems[i], algo)
#     end
#     EE = results[1:L+1]
#     TMI = results[L+2] + results[L+3] + results[L+4] - results[L+5] - results[L+6] - results[L+7] + results[L+8]
#     return EE, TMI
# end