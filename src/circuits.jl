abstract type Circuit end

struct YaoKivelsonCircuit <: Circuit
    size::Int
    nqubits::Int
    params::Vector{Any}
    equilibration_steps::Int64
    samples::Int64
    subsamples::Int64
    subsample_interval::Int64
    function YaoKivelsonCircuit(size, params; equilibration_steps=3*size, samples=100, subsamples=100, subsample_interval=1)
        return new(size, 3*size^2, params, equilibration_steps, samples, subsamples, subsample_interval)
    end
end