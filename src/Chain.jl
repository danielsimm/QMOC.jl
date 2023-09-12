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
    XX_projectors::Vector{PauliOperator}
    YY_projectors::Vector{PauliOperator}
    ZZ_projectors::Vector{PauliOperator}
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
    XX_projectors::Vector{PauliOperator}
    XY_projectors::Vector{PauliOperator}
    XZ_projectors::Vector{PauliOperator}
    YX_projectors::Vector{PauliOperator}
    YY_projectors::Vector{PauliOperator}
    YZ_projectors::Vector{PauliOperator}
    ZX_projectors::Vector{PauliOperator}
    ZY_projectors::Vector{PauliOperator}
    ZZ_projectors::Vector{PauliOperator}
end

include("ChainOperators.jl")

function _PPChainTrajectory(size::Int, nqubits::Int, name::String, params, checkpoints::Bool, verbosity::Symbol, index::Int64, thermalization_steps::Int64, measurement_steps::Int64, number_of_measurements::Int64) ::PPChainTrajectory
    XX_projectors = _chain_XX(nqubits)
    YY_projectors = _chain_YY(nqubits)
    ZZ_projectors = _chain_ZZ(nqubits)
    return PPChainTrajectory(size, nqubits, name, params, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements, XX_projectors, YY_projectors, ZZ_projectors)
end

function _PQChainTrajectory(size::Int, nqubits::Int, name::String, params, checkpoints::Bool, verbosity::Symbol, index::Int64, thermalization_steps::Int64, measurement_steps::Int64, number_of_measurements::Int64) ::PQChainTrajectory
    XX_projectors = _chain_XX(nqubits)
    XY_projectors = _chain_XY(nqubits)
    XZ_projectors = _chain_XZ(nqubits)
    YX_projectors = _chain_YX(nqubits)
    YY_projectors = _chain_YY(nqubits)
    YZ_projectors = _chain_YZ(nqubits)
    ZX_projectors = _chain_ZX(nqubits)
    ZY_projectors = _chain_ZY(nqubits)
    ZZ_projectors = _chain_ZZ(nqubits)
    return PQChainTrajectory(size, nqubits, name, params, checkpoints, verbosity, index, thermalization_steps, measurement_steps, number_of_measurements, XX_projectors, XY_projectors, XZ_projectors, YX_projectors, YY_projectors, YZ_projectors, ZX_projectors, ZY_projectors, ZZ_projectors)
end


function initialise(trajectory::ChainTrajectory) ::QuantumClifford.MixedDestabilizer
    return one(MixedDestabilizer, trajectory.size)
end

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::ChainTrajectory) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    for subtime in 1:trajectory.nqubits
        project!(state, random_operator(trajectory), phases=false, keep_result=false)
    end
end


### Observables

function entropy(state::QuantumClifford.MixedDestabilizer, trajectory::ChainTrajectory; algo=Val(:rref))
    N = trajectory.nqubits
    EE = zeros(33)
    subsystems = Int.(collect(range(0, N, 33)))
    for i in 2:32
        EE[i] = QuantumClifford.entanglement_entropy(state, 1:subsystems[i], algo)
    end
    return EE
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