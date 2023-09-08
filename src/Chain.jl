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

struct PPChainTrajectoryFast <: ChainTrajectory
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


function initialise(trajectory::ChainTrajectory)
    return one(MixedDestabilizer, trajectory.size)
end

function random_operator(trajectory::PPChainTrajectory)
    numberOfQubits = trajectory.nqubits
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]

    firstsite = rand(1:numberOfQubits)
    secondsite = mod1(firstsite+1, numberOfQubits)
    
    X_arr = falses(numberOfQubits)
    Z_arr = falses(numberOfQubits)

    probability = rand()

    if probability < px
        X_arr[firstsite] = true
        X_arr[secondsite] = true
    elseif probability < (px + py)
        X_arr[firstsite] = true
        X_arr[secondsite] = true
        Z_arr[firstsite] = true
        Z_arr[secondsite] = true
    else
        Z_arr[firstsite] = true
        Z_arr[secondsite] = true
    end

    return PauliOperator(0x00, X_arr,Z_arr)
end
    
function random_operator(trajectory::PQChainTrajectory)
    numberOfQubits = trajectory.nqubits
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]

    site1 = rand(1:numberOfQubits)
    site2 = mod1(firstsite+1, numberOfQubits)
    
    X_arr = falses(numberOfQubits)
    Z_arr = falses(numberOfQubits)

    probability1 = rand()
    if probability1 < px
        X_arr[site1] = true
    elseif probability1 < (px + py)
        X_arr[site1] = true
        Z_arr[site1] = true
    else
        Z_arr[site1] = true
    end
    probability2 = rand()
    if probability2 < px
        X_arr[site2] = true
    elseif probability2 < (px + py)
        X_arr[site2] = true
        Z_arr[site2] = true
    else
        Z_arr[site2] = true
    end

    return PauliOperator(0x00, X_arr,Z_arr)
end

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::ChainTrajectory) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    for subtime in 1:trajectory.nqubits
        project!(state, random_operator(trajectory), phases=false, keep_result=false)
    end
end

function circuit!(state::QuantumClifford.AbstractStabilizer, trajectory::PPChainTrajectoryFast) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]
    for subtime in 1:trajectory.nqubits
        probability = rand()
        if probability < px
            randomXX!(state, trajectory)
        elseif probability < (px + py)
            randomYY!(state, trajectory)
        else
            randomZZ!(state, trajectory)
        end
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