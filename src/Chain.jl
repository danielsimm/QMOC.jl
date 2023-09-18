abstract type ChainTrajectory <: Trajectory end

struct PPChainTrajectory <: ChainTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

struct PQChainTrajectory <: ChainTrajectory
    size::Int
    nqubits::Int
    name::String
    params::Vector{Real}
    checkpoints::Bool
    verbosity::Symbol
    index::Int64
    thermalization_steps::Int64
    measurement_steps::Int64
    number_of_measurements::Int64
end

include("ChainOperators.jl")


function initialise(trajectory::ChainTrajectory) ::QuantumClifford.MixedDestabilizer
    return one(MixedDestabilizer, trajectory.size)
end

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::ChainTrajectory, operators) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    for subtime in 1:trajectory.nqubits
        project!(state, random_operator(trajectory, operators), phases=false, keep_result=false)
    end
end


### Observables

function entropy(state::QuantumClifford.MixedDestabilizer, trajectory::ChainTrajectory; algo=Val(:rref))
    N = trajectory.nqubits
    EE = zeros(33)
    subsystems = round.(Int, collect(range(0, N, 33)))
    for i in 2:32
        EE[i] = QuantumClifford.entanglement_entropy(state, 1:subsystems[i], algo)
    end
    return EE
end

function subsystem_labels(trajectory::ChainTrajectory)
    N = trajectory.nqubits
    return round.(Int, collect(range(0, N, 33)))
end

function tmi(state::QuantumClifford.MixedDestabilizer, trajectory::ChainTrajectory; algo=Val(:rref)) # no geometry
    N = trajectory.nqubits
    A = 1:Int(N/4)
    B = Int(N/4)+1:Int(N/2)
    C = Int(N/2)+1:Int(3N/4)
    
    if mod(N, 4) != 0
        @info "L must be a multiple of 4, but is $(N). Tripartite mutual information is ill-defined."
    end

    SA = entanglement_entropy(state, A, algo)
    SB = entanglement_entropy(state, B, algo)
    SC = entanglement_entropy(state, C, algo)
    SAB = entanglement_entropy(state, union(A,B), algo)
    SBC = entanglement_entropy(state, union(B,C), algo)
    SAC = entanglement_entropy(state, union(A,C), algo)
    SABC = entanglement_entropy(state, union(A, B, C), algo)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end