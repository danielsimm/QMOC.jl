using QMOC
using DelimitedFiles

if !isdir("kitaev_phaseboundary")
    mkdir("kitaev_phaseboundary")
end
cd("kitaev_phaseboundary")
starting_points = QMOC.parameter_line(:px, [1//2, 1//2, 0], 15)
endpoint = [1//3, 1//3, 1//3]

lines = [QMOC.parameter_line(starting_points[i], endpoint, 15) for i in eachindex(starting_points)]
Ls = [12, 20, 28, 32, 40]

for i in eachindex(Ls)
    L = Ls[i]
    if !isdir("L=$(L)")
        mkdir("L=$(L)")
    end
    for j in eachindex(lines)
        println("L=$(L), line=$(j)")
        parameter_set = lines[j]
        data = zeros(length(parameter_set), 4)
        for p in eachindex(parameter_set)
            params = parameter_set[p]
            data[p, 1:3] = params
            c = QMOC.KitaevCircuit(L, params)
            operators = QMOC.get_operators(c)
            n_states = 18
            n_samples = 200
            this_I3 = zeros(n_states)
            Threads.@threads for st in 1:n_states
                state = QMOC.thermal_state(c)
                for _ in 1:n_samples
                    QMOC.apply!(state, c, operators)
                    this_I3[st] += QMOC.tmi(state, c)
                end
                this_I3[st] /= n_samples
            end
            data[p, 4] = mean(this_I3)
        end
        writedlm("L=$(L)/line_$(j).txt", data)
    end
end
