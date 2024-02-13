abstract type Circuit end

struct YaoKivelsonCircuit <: Circuit
    size::Int
    nqubits::Int
    operators::Vector{QuantumClifford.PauliOperator}
    weights::Vector{Float64}
    observer::Function
    if length(operators) != length(weights)
        throw(ArgumentError("The number of operators and weights must be the same"))
    end
    function YaoKivelsonCircuit(size, params; equilibration_steps=3*size, samples=100, subsamples=100, subsample_interval=1)
        return new(size, 3*size^2, params, equilibration_steps, samples, subsamples, subsample_interval)
    end
end

function apply!(state, operators, circuit)