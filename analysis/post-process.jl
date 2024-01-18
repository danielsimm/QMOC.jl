using QMOC
using DelimitedFiles


function write_out(name)
    if isfile("../LoopPaperCliffordResults/data/$(name).txt")
        return nothing
    end
    I3 = Float64.(vec(QMOC.evaluate(name, :I3)[2]))
    I3_err = Float64.(vec(QMOC.evaluate(name, :ΔI3)[2])) ./ sqrt(QMOC.loadMetadata(name).num_trajectrories)
    params = QMOC.loadMetadata(name).parameter_set
    Js = [params[i][1] for i in eachindex(params)]
    Ks = [params[i][2] for i in eachindex(params)]
    writedlm("../LoopPaperCliffordResults/data/$(name).txt", hcat(I3, I3_err, Js, Ks))
    return nothing
end

sims = [
    "nonorientable_scaling_24",
    "nonorientable_scaling_32",
    "orientable_scaling_16",
    "orientable_scaling_24",
    "orientable_scaling_32",
    "orientable_scaling_40",
    "orientable_scaling_48",
]

write_out.(sims)