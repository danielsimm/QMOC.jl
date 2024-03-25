using DelimitedFiles
using CairoMakie

L = 24
function K_range(L, nu, Kc, ΔK)
    K_min = max(0, Kc - ΔK/(L^(1/nu)))
    K_max = Kc + ΔK/(L^(1/nu))
    return LinRange(K_min, K_max, 30)
end
Kc_orientable = 0.631
nu_orientable = 1.0
Kc_nonorientable = 0.654
nu_nonorientable = 0.94
ΔK = 4
Ks_orientable = K_range(L, nu_orientable, Kc_orientable, ΔK)
Ks_nonorientable = K_range(L, nu_nonorientable, Kc_nonorientable, ΔK)
parameter_set_orientable = [[1-K, K] for K in Ks_orientable]
parameter_set_nonorientable = [[1-K, K] for K in Ks_nonorientable]

avg = readdlm("tmis_YaoKivelsonOrientable_24.txt")
err = readdlm("errors_YaoKivelsonOrientable_24.txt")
x = (0:1:length(avg[1, :])-1) ./ (L)
ps = [K for K in Ks_orientable]
with_theme(theme_latexfonts()) do 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1])
    ax.title = "Equilibration times"
    for i in eachindex(ps)
        scatterlines!(ax, x, avg[i, :], label="p = $(ps[i])")
        errorbars!(ax, x, avg[i, :], err[i, :])
    end
    #axislegend(ax, position = :rc)
    fig
end