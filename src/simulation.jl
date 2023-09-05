include("LatticeCircuits.jl")
Ls = [12]
parameter_set = parameter_full(10)

#usage: simulate(L, model, 1:max_trajectory, parameter_set, index, name, checkpoints=false, verbose=false)
for L in Ls
    for index in eachindex(parameter_set)
        simulate(L, :YaoKivelsonXYZ, 16, parameter_set, index, "full", false, true)
    end
end