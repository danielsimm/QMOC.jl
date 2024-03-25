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

include("KitaevKekuleOperators.jl")


function get_operators(trajectory::KitaevTrajectory)
    matrix = Matrix{PauliOperator}(undef, 3, trajectory.size^2)
    matrix[1, :] .= _HC_XX_operators(trajectory.size)
    matrix[2, :] .= _HC_YY_operators(trajectory.size)
    matrix[3, :] .= _HC_ZZ_operators(trajectory.size)
    return matrix
end

function get_operators(trajectory::KekuleTrajectory)
    matrix = Matrix{PauliOperator}(undef, 3, trajectory.size^2)
    matrix[1, :] .= _HC_red_operators(trajectory.size)
    matrix[2, :] .= _HC_green_operators(trajectory.size)
    matrix[3, :] .= _HC_blue_operators(trajectory.size)
    return matrix
end
    
function initialise(trajectory::HoneycombTrajectory)
    L = trajectory.size
    stabs = [_HC_ZZ_operators(L)[1:L^2-1]..., _HC_WilsonPlaquette_operators(L)..., _HC_WilsonLoops(L)...]
    stabiliser = Stabilizer(stabs)
    if QuantumClifford.trusted_rank(stabiliser) != 2*L^2
        @warn "Initial state is not pure."
    end
    return MixedDestabilizer(stabiliser)
end


### Dynamics ###

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::KitaevTrajectory, operators)
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]
    L = trajectory.size
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < px
            project!(state, operators[1, rand(1:L^2)], keep_result=false, phases=false)
        elseif p < px + py
            project!(state, operators[2, rand(1:L^2)], keep_result=false, phases=false)
        elseif p < px + py + pz
            project!(state, operators[3, rand(1:L^2)], keep_result=false, phases=false)
        end
    end
    return nothing
end

function circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::KekuleTrajectory, operators)
    pr = trajectory.params[1]
    pg = trajectory.params[2]
    pb = trajectory.params[3]
    L = trajectory.size
    for subtime in 1:trajectory.nqubits
        p = rand()
        if p < pr
            project!(state, operators[1, rand(1:L^2)], keep_result=false, phases=false)
        elseif p < pr + pg
            project!(state, operators[2, rand(1:L^2)], keep_result=false, phases=false)
        elseif p < pr + pg + pb
            project!(state, operators[3, rand(1:L^2)], keep_result=false, phases=false)
        end
    end
    return nothing
end

function apply!(stabilizer, c, operators)
    p1 = c.params[1]
    p2 = c.params[2] + p1
    N = c.size^2
    for subtime in 1:c.nqubits
        p = rand()
        if p < p1
            project!(stabilizer, operators[1, rand(1:N)], keep_result=false, phases=false)
        elseif p < p2
            project!(stabilizer, operators[2, rand(1:N)], keep_result=false, phases=false)
        else
            project!(stabilizer, operators[3, rand(1:N)], keep_result=false, phases=false)
        end
    end
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

# function pure_honeycomb_entanglement(state::QuantumClifford.MixedDestabilizer, subsystem_qubits, num_qubits, L)
#     stab = stabilizerview(state)
#     N = num_qubits
#     n_crossings = 0
#     non_subsystem_qubits = setdiff(1:N, subsystem_qubits)
#     for i in 1:N
#         non_idents = xbit(stab[i]) .|| zbit(stab[i])
#         non_idents_subA = non_idents[subsystem_qubits]
#         non_idents_subB = non_idents[non_subsystem_qubits]
#         if sum(non_idents_subA) > 0 && sum(non_idents_subB) > 0
#             n_crossings += 1
#         end
#     end
#     return L/2 - 1 + n_crossings/2
# end

    

# function new_entropy(state::QuantumClifford.MixedDestabilizer, trajectory::HoneycombTrajectory)
#     L = trajectory.size
#     EE = zeros(L+1)
#     for i in 1:L
#         EE[i+1] = pure_honeycomb_entanglement(state, HC_subsystem(L, 1:i), trajectory.nqubits, L)
#     end
#     return EE
# end

function subsystem_labels(trajectory::HoneycombTrajectory)
    L = trajectory.size
    subsystems = zeros(L+1)
    for i in 1:L
        subsystems[i+1] = i
    end
    return subsystems
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