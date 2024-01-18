using QMOC
using DelimitedFiles


function write_out(name)
    if isfile("../ProjectiveMeasurementCircuits.jl//data/$(name).txt")
        return nothing
    end
    I3 = Float64.(vec(QMOC.evaluate(name, :I3)[2]))
    I3_err = Float64.(vec(QMOC.evaluate(name, :ΔI3)[2])) ./ sqrt(QMOC.loadMetadata(name).num_trajectrories)
    params = QMOC.loadMetadata(name).parameter_set
    Js = [params[i][1] for i in eachindex(params)]
    Ks = [params[i][2] for i in eachindex(params)]
    writedlm("../ProjectiveMeasurementCircuits.jl/data/$(name).txt", hcat(Js, Ks, I3, I3_err))
    return nothing
end

sims = [
    "nonorientable_scaling_40",
    "orientable_scaling_40",
]

write_out.(sims)