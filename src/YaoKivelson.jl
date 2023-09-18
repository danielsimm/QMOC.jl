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

include("YaoKivelsonOperators.jl")

function get_operators(trajectory::YaoKivelsonXYZTrajectory)
    matrix = Matrix{PauliOperator}(undef, 3, 3*trajectory.size^2)
    matrix[1, :] = _DHC_XX_operators(trajectory.size)
    matrix[2, :] = _DHC_YY_operators(trajectory.size)
    matrix[3, :] = _DHC_ZZ_operators(trajectory.size)
    return matrix
end

function get_operators(trajectory::YaoKivelsonJJTrajectory)
    matrix = Matrix{PauliOperator}(undef, 2, 3*trajectory.size^2)
    matrix[1, :] = _DHC_J_operators(trajectory.size)
    matrix[2, :] = _DHC_Jprime_operators(trajectory.size)
    return matrix
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


### Dynamics ###

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonXYZTrajectory, operators)
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]
    L = trajectory.size
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < px
            project!(state, operators[1, rand(1:3*L^2)], keep_result=false, phases=false)
        elseif p < px + py
            project!(state, operators[2, rand(1:3*L^2)], keep_result=false, phases=false)
        elseif p < px + py + pz
            project!(state, operators[3, rand(1:3*L^2)], keep_result=false, phases=false)
        end
    end
    return nothing
end

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonJJTrajectory, operators)
    pJ = trajectory.params[1]
    pJprime = trajectory.params[2]
    L = trajectory.size
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < pJ
            project!(state, operators[1, rand(1:3*L^2)], keep_result=false, phases=false)
        elseif p < pJ + pJprime
            project!(state, operators[2, rand(1:3*L^2)], keep_result=false, phases=false)
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

function subsystem_labels(trajectory::DecoratedHoneycombTrajectory)
    L = trajectory.size
    subsystems = zeros(L+1)
    for i in 1:L
        subsystems[i+1] = i
    end
    return subsystems
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