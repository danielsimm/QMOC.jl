abstract type DecoratedHoneycombTrajectory <: Trajectory end

abstract type YaoKivelsonJKTrajectory <: DecoratedHoneycombTrajectory end
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

struct YaoKivelsonOrientableTrajectory <: YaoKivelsonJKTrajectory
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

struct YaoKivelsonNonorientableTrajectory <: YaoKivelsonJKTrajectory
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

export YaoKivelsonXYZTrajectory, YaoKivelsonNonorientableTrajectory, YaoKivelsonOrientableTrajectory, DecoratedHoneycombTrajectory

### Operators ###

include("YaoKivelsonOperators.jl")

function get_operators(trajectory::YaoKivelsonXYZTrajectory)
    matrix = Matrix{PauliOperator}(undef, 3, 3*trajectory.size^2)
    matrix[1, :] = _DHC_XX_operators(trajectory.size)
    matrix[2, :] = _DHC_YY_operators(trajectory.size)
    matrix[3, :] = _DHC_ZZ_operators(trajectory.size)
    return matrix
end

function get_operators(trajectory::YaoKivelsonOrientableTrajectory)
    L = trajectory.size
    vector = Vector{PauliOperator}(undef, (3+4)*L^2)
    vector[1:4*L^2] = _DHC_J_operators_orientable(L)
    vector[4*L^2+1:end] = _DHC_K_operators(L)
    return vector
end

function get_operators(trajectory::YaoKivelsonNonorientableTrajectory)
    L = trajectory.size
    vector = Vector{PauliOperator}(undef, (3+6)*L^2)
    vector[1:6*L^2] = _DHC_J_operators_nonorientable(L)
    vector[6*L^2+1:end] = _DHC_K_operators(L)
    return vector
end

# function initialise_mixed(trajectory::DecoratedHoneycombTrajectory; basis=:Z) :: MixedDestabilizer
#     L = trajectory.size
#     if basis == :Z
#         bilinears = _DHC_ZZ_operators(L)
#     elseif basis == :X
#         bilinears = _DHC_XX_operators(L)
#     elseif basis == :Y
#         bilinears = _DHC_YY_operators(L)
#     end
#     stabs = [bilinears..., _DHC_largeloop_operators(L)..., _DHC_smallloop_operators(L)..., _DHC_wilsonline_operators(L)...]
#     state = MixedDestabilizer(Stabilizer(stabs))
#     if (QuantumClifford.trusted_rank(state) != trajectory.nqubits) && (trajectory.verbosity == :debug)
#         @warn "Initial state is not pure."
#     end
#     return state
# end

function initialise(trajectory::DecoratedHoneycombTrajectory)
    L = trajectory.size
    stab = QuantumClifford.MixedDestabilizer(QuantumClifford.Stabilizer(one(QuantumClifford.Tableau, 6*L^2; basis=:X)))
    largeloops = _DHC_largeloop_operators(L)
    smallloops = _DHC_smallloop_operators(L)
    wilsonlines = _DHC_wilsonline_operators(L)
    for op in [largeloops..., smallloops..., wilsonlines...]
        QuantumClifford.project!(stab, op, phases=false)
    end
    @assert QuantumClifford.trusted_rank(stab) == 6*L^2
    return stab
end

### Dynamics ###

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonXYZTrajectory, operators::Vector{PauliOperator}) ::Nothing
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

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::YaoKivelsonOrientableTrajectory, operators::Vector{PauliOperator}) ::Nothing
    J = trajectory.params[1] # probability of triangular measurement, orientable: silence triangular Z bond, non-orientable: all bonds
    # K = trajectory.params[2] # probability of hexagonal/Kitaev-style measurement
    for subtime in 1:trajectory.nqubits
        if rand() < J
            project!(state, operators[rand(1:4*trajectory.size^2)], keep_result=false, phases=false)
        else
            project!(state, operators[4*trajectory.size^2+rand(1:3*trajectory.size^2)], keep_result=false, phases=false)
        end
    end
    return nothing
end

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::YaoKivelsonNonorientableTrajectory, operators::Vector{PauliOperator}) ::Nothing
    J = trajectory.params[1] # probability of triangular measurement, orientable: silence triangular Z bond, non-orientable: all bonds
    # K = trajectory.params[2] # probability of hexagonal/Kitaev-style measurement
    for subtime in 1:trajectory.nqubits
        if rand() < J
            project!(state, operators[rand(1:6*trajectory.size^2)], keep_result=false, phases=false)
        else
            project!(state, operators[6*trajectory.size^2+rand(1:3*trajectory.size^2)], keep_result=false, phases=false)
        end
    end
    return nothing
end


### Observables ###

"""
    DHC_snake(L)

Return a vector of sites in the order of the snake pattern.
"""
function DHC_snake(L)
    sites = zeros(Int, 6*L^2)
    rows = [((i-1)*6*L)+1:i*6*L for i in 1:L]
    for r in eachindex(rows)
        if iseven(r)
            sites[rows[r]] = rows[r][end:-1:1]
        else
            sites[rows[r]] = rows[r]
        end
    end
    return sites
end

function new_entropy(state, trajectory::DecoratedHoneycombTrajectory)
    stab = copy(stabilizerview(state))
    L = trajectory.size
    QuantumClifford.permute!(stab, DHC_snake(L))
end


function entropy(state::QuantumClifford.AbstractStabilizer, trajectory::DecoratedHoneycombTrajectory)
    algo=Val(:rref)
    L = trajectory.size
    EE = zeros(L+1)
    # for i in 1:L
    #     EE[i+1] = entanglement_entropy(state, DHC_subsystem(L, 1:i), algo)
    # end
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

function tmi(state::QuantumClifford.AbstractStabilizer, trajectory::DecoratedHoneycombTrajectory)
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

function YaoKivelsonTest(L, mode, J)
    return trajectory(
        mode,
        L,
        6*L^2,
        "YK_$(mode)_test",
        [J, 1-J],
        false,
        :debug,
        1,
        3*L,
        1,
        1
    )
end

function YaoKivelsonNonorientableTest(L)
    return trajectory(
        :YaoKivelsonNonorientable,
        L,
        6*L^2,
        "YK_nonorientable_test",
        [1-0.654, 0.654],
        false,
        :debug,
        1,
        3*L,
        1,
        1
    )
end