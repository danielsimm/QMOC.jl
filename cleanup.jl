using JLD2
using DelimitedFiles
folder = "cluster/Kitaev28"
files = readdir(folder)
to_read = filter(contains("idx"), files)
idx_array = sort(unique(parse.(Int64, [split(split(string, "idx")[2], "_")[1] for string in to_read])))

Threads.@threads for idx in idx_array
    this_files = filter(contains("idx$(idx)_"), to_read)
    out = zeros(length(this_files))
    for i in eachindex(this_files)
        file = this_files[i]
        out[i] = jldopen("$(folder)/$(file)")["I3"]
        rm("$(folder)/$(file)")
    end
    writedlm("$(folder)/idx$(idx).txt", out)
end

folder = "cluster/Kitaev32"
files = readdir(folder)
to_read = filter(contains("idx"), files)
idx_array = sort(unique(parse.(Int64, [split(split(string, "idx")[2], "_")[1] for string in to_read])))

Threads.@threads for idx in idx_array
    this_files = filter(contains("idx$(idx)_"), to_read)
    out = zeros(length(this_files))
    for i in eachindex(this_files)
        file = this_files[i]
        out[i] = jldopen("$(folder)/$(file)")["I3"]
        rm("$(folder)/$(file)")
    end
    writedlm("$(folder)/idx$(idx).txt", out)
end