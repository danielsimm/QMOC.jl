using QMOC
using CairoMakie

params, I3_NO_12 =evaluate("YaoKivelsonNonorientableFine_12", :I3)
I3_NO_24 = evaluate("YaoKivelsonNonorientableFine_24", :I3)[2]
I3_NO_36 = evaluate("YaoKivelsonNonorientableFine_36", :I3)[2]
I3_NO_12_err = evaluate("YaoKivelsonNonorientableFine_12", :ΔI3)[2]
I3_NO_24_err = evaluate("YaoKivelsonNonorientableFine_24", :ΔI3)[2]
I3_NO_36_err = evaluate("YaoKivelsonNonorientableFine_36", :ΔI3)[2]
I3_NO_12 = Float64.(I3_NO_12)
I3_NO_24 = Float64.(I3_NO_24)
I3_NO_36 = Float64.(I3_NO_36)
I3_NO_12_err = Float64.(I3_NO_12_err)
I3_NO_24_err = Float64.(I3_NO_24_err)
I3_NO_36_err = Float64.(I3_NO_36_err)
blue = colorant"#175EC7"
red = colorant"#C81E5B"
colormap = ColorScheme(range(red, blue))
c1 = get(colormap, 0.0)
c2 = get(colormap, 0.6)
c3 = get(colormap, 1.0)
with_theme(merge(theme_black(), theme_latexfonts())) do
J = [params[i][1] for i in eachindex(params)]
K = [params[i][2] for i in eachindex(params)]
x = K./J
fig = Figure(resolution = (700, 400), backgroundcolor = (:black, 0.0), textcolor=:white, axiscolor=:white, color=:white)
ax = Axis(fig[1, 1], xlabel = L"\frac{K}{J}", ylabel = L"I(A:B:C)", backgroundcolor=(:black,0.0))
lw = 1.8
err_lw = 1.1
ms = 8
p3 = scatterlines!(ax, x[2:end], I3_NO_12[2:end], linewidth = lw, markersize = ms, label = "L = 12", color=c3)
scatter!(ax, x[1], I3_NO_12[1], markersize = ms, color=c3)
errorbars!(ax, x, I3_NO_12, I3_NO_12_err, color=c3, linewidth=err_lw)

p2 = scatterlines!(ax, x[2:end], I3_NO_24[2:end], linewidth = lw, markersize = ms, label = "L = 24", color=c2)
scatter!(ax, x[1], I3_NO_24[1], markersize = ms, color=c2)
errorbars!(ax, x, I3_NO_24, I3_NO_24_err, color=c2, linewidth=err_lw)

p1 = scatterlines!(ax, x[2:end], I3_NO_36[2:end], linewidth = lw, markersize = ms, label = "L = 36", color=c1)
scatter!(ax, x[1], I3_NO_36[1], markersize = ms, color=c1)
errorbars!(ax, x, I3_NO_36, I3_NO_36_err, color=c1, linewidth=err_lw)

Legend(fig[1,2], [p1, p2, p3], ["L=36", "L=24", "L=12"], bgcolor=(:black, 0.0))
save("YaoKivelsonNonorientable_coarse.pdf", fig)
end