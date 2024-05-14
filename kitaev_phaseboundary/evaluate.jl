using QMOC
using ScalingCollapse
Ls = [12, 20, 28, 32]
cd("kitaev_phaseboundary")
startpoints = parameter_line(:px, [1//2, 1//2, 0], 15)
critical_distance = zeros(length(startpoints))

for i in eachindex(startpoints)
    line = parameter_line(startpoints[i], [1//3, 1//3, 1//3], 15)
    distances = [QMOC.parameter_distance.(line) for _ in Ls]
    I3s =[zeros(length(line)) for _ in Ls]
    for j in eachindex(Ls)
        data = readdlm("L=$(Ls[j])/line_$(i).txt")
        I3s[j] = data[:, 4]
    end
    sp = ScalingProblem(distances, I3s, Ls, sf=ScalingFunction(:x),)
    critical_distance[i] = sp.optimal_ps[1]
end
critical_distance
critical_params = zeros(length(startpoints), 3)
for i in eachindex(startpoints)
    line = parameter_line(startpoints[i], [1//3, 1//3, 1//3], 1000)
    distances = QMOC.parameter_distance.(line)
    crit_ind = argmin(abs.(distances .- critical_distance[i]))
    critical_params[i, :] = line[crit_ind]
end

critical_params

cartesian = [QMOC.parametric_to_cartesian(critical_params[i, :]) for i in 1:size(critical_params, 1)]

using CairoMakie

begin
    fig = Figure()
    ax = Axis(fig[1,1])
    i = 1
    for L in Ls
        data = readdlm("L=$(L)/line_$(i).txt")
        x = [QMOC.parameter_distance(data[j, 1:3]) for j in 1:size(data, 1)]
        y = data[:, 4]
        scatterlines!(ax, x, y, label="L=$(L)")
    end
    fig
end

begin
    fig = Figure()
    ax = Axis(fig[1,1])
    for i in 1:size(critical_params, 1)
        scatter!(ax, [cartesian[i][1]], [cartesian[i][2]], markersize=10)
    end
    fig
end

writedlm("critical_params.txt", critical_params)