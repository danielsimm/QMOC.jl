using DelimitedFiles
using CairoMakie

avg = readdlm("equilibration_times_avg.txt")
err = readdlm("equilibration_times_err.txt")
x = (0:1:length(avg[1, :])-1) ./ (12)
ps = [0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0, "p_c"]
with_theme(theme_latexfonts()) do 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1])
    ax.title = "Equilibration times"
    for i in 1:8
        scatterlines!(ax, x, avg[i, :], label="p = $(ps[i])")
        errorbars!(ax, x, avg[i, :], err[i, :])
    end
    axislegend(ax, position = :rc)
    fig
end