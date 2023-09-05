using QuantumClifford
using Statistics
using Dates
using LinearAlgebra

struct Setup
    mode::Symbol
    size::Int
    name::String
    params::Vector{Vector{Real}}
    trajectories::UnitRange{Int64}
    checkpoints::Bool
    verbosity::Symbol
end

function simulate(setup::Setup)
end

function _id(setup::Setup, index::Int64, trajectory::Int64)
    param = setup.params[index]
    param_string = ""
    for p in param
        param_string *= "$(p)_"
    end
    return "$(setup.mode)_$(setup.name)_$(setup.size)_$(trajectory)"*param_string
end

# function _update_checkpoint(state, id, time)
#     filename = "data/checkpoints/$(hash(id)).jld2" 
#     if isfile(filename)
#         if time < jldopen(filename, "r") do file
#             file["time"]
#         end
#             jldopen(filename, "a+") do file
#                 file["state"] = state
#                 file["time"] = time
#             end
#         end
#         # load
#         # return state
#     # else
#         # (@async) save checkpoint (state, id, time) -> overwrite
#         # return nothing
# end

# function _filehandling(setup::Setup)
#     if !(isempty("data/$(setup.mode)/$(setup.name)/L=$(setup.)")

#     # if !ispath("data/$(setup.type)")
#     #     mkpath("data/$(setup.type)")
#     # end
#     # if !ispath("data/$(setup.type)/$(setup.size)")
#     #     mkpath("data/$(setup.type)/$(setup.size)")
#     # end
#     # if !ispath("data/$(setup.type)/$(setup.size)/$(setup.trajectories)")
#     #     mkpath("data/$(setup.type)/$(setup.size)/$(setup.trajectories)")
#     # end
# end
function simulate(L, mode, trajectories, parameter_set, index, name, checkpoints=false, verbose=false)
    # set BLAS threads to 1 to avoid oversubscription
    BLAS.set_num_threads(1)
    
    parameters = parameter_set[index]
    if !ispath("data/$(mode)/$(name)/L=$(L)")
        mkpath("data/$(mode)/$(name)/L=$(L)")
    end
    if !ispath("data/checkpoints")
        mkpath("data/checkpoints")
    end
    if isfile("data/$(mode)/$(name)/L=$(L)/param_$(index).jld2")
        @info "Skipping $(mode) | $(name) | L=$(L) | Parameter Set $(index), already exists."
        return
    else
        if mode == :ChainPP
            simulate_chain(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        elseif mode == :ChainPQ
            simulate_chain(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        elseif mode == :Kekule
            simulate_honeycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        elseif mode == :Kitaev
            simulate_honeycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        elseif mode == :YaoKivelsonJJ
            simulate_decoratedhoneycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        elseif mode == :YaoKivelsonXYZ
            simulate_decoratedhoneycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
        else
            error("Mode $(mode) not implemented.")
        end
    end
end

function simulate_decoratedhoneycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
    if typeof(trajectories) == Int64
        trajectories = collect(1:trajectories)
        verbose ? (@info "Converting to trajectory range $(trajectories).") : nothing
    end
    
    if isfile("data/$(mode)/init_L=$(L).jld2")
        init = jldopen("data/$(mode)/init_L=$(L).jld2") do file
            file["init"]
        end
    else
        init = DHC_initial_state(L)
        jldopen("data/$(mode)/init_L=$(L).jld2", "a+") do file
            file["init"] = init
        end
        verbose ? (@info "Saved initial state to file for L=$(L) decorated honeycomb lattice.") : nothing
    end

    @info "Starting $(mode) | $(name) | L=$(L) | parameter set $(index)."
    tick = now()

    entropy = zeros(length(trajectories), L+1)
    TMI = zeros(length(trajectories))
    Threads.@threads for traj in trajectories
        id = "$(mode)_$(name)_$(L)_$(index)_$(traj)"
        entropy[traj, :], TMI[traj] = trajectory(id, copy(init), L, parameters, mode, checkpoints, verbose)
        
    end

    jldopen("data/$(mode)/$(name)/L=$(L)/param_$(index).jld2", "a+") do file
        file["entropy"] = mean(entropy, dims=1)
        file["TMI"] = mean(TMI)
        file["entropy_std"] = std(entropy, dims=1)
        file["TMI_std"] = std(TMI)
    end
    tock = now()
    @info "Finished $(mode) | $(name) | L=$(L) | parameter set $(index). Time elapsed: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS"))."

    # remove all checkpoints
    if checkpoints
        for traj in trajectories
            rm("data/checkpoints/$(mode)_$(name)_$(L)_$(index)_$(traj).jld2")
        end
    end
end

function simulate_honeycomb(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
    if typeof(trajectories) == Int64
        trajectories = collect(1:trajectories)
        verbose ? (@info "Converting to trajectory range $(trajectories).") : nothing
    end
    
    if isfile("data/$(mode)/init_L=$(L).jld2")
        init = jldopen("data/$(mode)/init_L=$(L).jld2") do file
            file["init"]
        end
    else
        init = HC_initial_state(L)
        jldopen("data/$(mode)/init_L=$(L).jld2", "a+") do file
            file["init"] = init
        end
        verbose ? (@info "Saved initial state to file for L=$(L) honeycomb lattice.") : nothing
    end

    @info "Starting $(mode) | $(name) | L=$(L) | parameter set $(index)."
    tick = now()

    entropy = zeros(length(trajectories), L+1)
    TMI = zeros(length(trajectories))
    Threads.@threads for traj in trajectories
        id = "$(mode)_$(name)_$(L)_$(index)_$(traj)"
        entropy[traj, :], TMI[traj] = trajectory(id, copy(init), L, parameters, mode, checkpoints, verbose)
        
    end

    jldopen("data/$(mode)/$(name)/L=$(L)/param_$(index).jld2", "a+") do file
        file["entropy"] = mean(entropy, dims=1)
        file["TMI"] = mean(TMI)
        file["entropy_std"] = std(entropy, dims=1)
        file["TMI_std"] = std(TMI)
    end
    tock = now()
    @info "Finished $(mode) | $(name) | L=$(L) | parameter set $(index). Time elapsed: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS"))."

    # remove all checkpoints
    if checkpoints
        for traj in trajectories
            rm("data/checkpoints/$(mode)_$(name)_$(L)_$(index)_$(traj).jld2")
        end
    end
end

function simulate_chain(L, mode, trajectories, parameters, index, name, checkpoints, verbose)
    init = one(MixedDestabilizer, L, L) 

    if typeof(trajectories) == Int64
        trajectories = 1:trajectories
        verbose ? (@info "Converting to trajectory range $(trajectories).") : nothing
    end
    

    @info "Starting $(mode) | $(name) | L=$(L) | parameter set $(index)."
    tick = now()

    entropy = zeros(length(trajectories), 33)
    TMI = zeros(length(trajectories))
    Threads.@threads for traj in trajectories
        id = "$(mode)_$(name)_$(L)_$(index)_$(trajectory)"
        entropy[traj, :], TMI[traj] = trajectory(id, copy(init), L, parameters, mode, checkpoints, verbose)
        
    end

    jldopen("data/$(mode)/$(name)/L=$(L)/param_$(index).jld2", "a+") do file
        file["entropy"] = mean(entropy, dims=1)
        file["TMI"] = mean(TMI)
        file["entropy_std"] = std(entropy, dims=1)
        file["TMI_std"] = std(TMI)
    end
    tock = now()
    @info "Finished $(mode) | $(name) | L=$(L) | parameter set $(index). Time elapsed: $(Dates.format(convert(DateTime, tock-tick), "HH:MM:SS"))."
end