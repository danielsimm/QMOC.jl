using QMOC
using CairoMakie
using ColorSchemes
using Scaling

begin
#loadMetadata("YaoKivelsonOrientableDetail_12").parameter_set[1:21] == loadMetadata("YaoKivelsonOrientableCritical_24").parameter_set
I3_12 = evaluate("YaoKivelsonOrientableDetail_12", :I3)[2][1:21]
params, I3_24 =evaluate("YaoKivelsonOrientableCritical_24", :I3)
I3_28 = evaluate("YaoKivelsonOrientableCritical_28", :I3)[2]
I3_32 = evaluate("YaoKivelsonOrientableCritical_32", :I3)[2]
I3_36 =evaluate("YaoKivelsonOrientableCritical_36", :I3)[2]
I3_40 = evaluate("YaoKivelsonOrientableCritical_40", :I3)[2]
I3_48 =evaluate("YaoKivelsonOrientableCriticalBetter_48", :I3)[2]
I3_60 = evaluate("YaoKivelsonOrientableCriticalBetter_60", :I3)[2]
# I3_72 = evaluate("YaoKivelsonOrientableCriticalFast_72", :I3)[2]
I3_12_err = evaluate("YaoKivelsonOrientableDetail_12", :ΔI3)[2][1:21]
I3_24_err = evaluate("YaoKivelsonOrientableCritical_24", :ΔI3)[2]
I3_28_err = evaluate("YaoKivelsonOrientableCritical_28", :ΔI3)[2]
I3_32_err = evaluate("YaoKivelsonOrientableCritical_32", :ΔI3)[2]
I3_36_err = evaluate("YaoKivelsonOrientableCritical_36", :ΔI3)[2]
I3_40_err = evaluate("YaoKivelsonOrientableCritical_40", :ΔI3)[2]  
I3_48_err = evaluate("YaoKivelsonOrientableCriticalBetter_48", :ΔI3)[2]
I3_60_err = evaluate("YaoKivelsonOrientableCriticalBetter_60", :ΔI3)[2]
# I3_72_err = evaluate("YaoKivelsonOrientableCriticalFast_72", :ΔI3)[2]

I3_12 = Float64.(I3_12)
I3_24 = Float64.(I3_24)
I3_28 = Float64.(I3_28)
I3_32 = Float64.(I3_32)
I3_36 = Float64.(I3_36)
I3_40 = Float64.(I3_40)
I3_48 = Float64.(I3_48)
I3_60 = Float64.(I3_60)
# I3_72 = Float64.(I3_72)
I3_12_err = Float64.(I3_12_err)./sqrt(200)
I3_24_err = Float64.(I3_24_err)./sqrt(200)
I3_28_err = Float64.(I3_28_err)./sqrt(200)
I3_32_err = Float64.(I3_32_err)./sqrt(200)
I3_36_err = Float64.(I3_36_err)./sqrt(200)
I3_40_err = Float64.(I3_40_err)./sqrt(100)
I3_48_err = Float64.(I3_48_err)./sqrt(100)
I3_60_err = Float64.(I3_60_err)./sqrt(100)
# I3_72_err = Float64.(I3_72_err)
end


begin
    x = Float64.([params[i][2] for i in eachindex(params)])
    data =[I3_12 I3_24 I3_28 I3_32 I3_36 I3_40 I3_48]
    err = [I3_12_err I3_24_err I3_28_err I3_32_err I3_36_err I3_40_err I3_48_err]
    for i in eachindex(err)
        if err[i] == 0
            err[i] = 1e-10
        end
    end
    perm = sortperm(x)
    x = x[perm]
    data = data[perm, :]
    err = err[perm, :]
end

sp = ScalingProblem(x, data[:, 1:end], [12, 24, 28, 32, 36, 40, 60], algorithm=:spline;
    sf=ScalingFunction(:x),
    dx=[-2.5,2.5],
    p_space = [0.5:0.1:0.7, 0.9:0.001:1.1],
    #verbose=true,
    error=true
)

p_space, residuals = Scaling.residual_landscapes(sp, p_space=[0.5:0.1:0.7, 0.9:0.001:1.1], algorithm=:spline)
# with_theme(merge(theme_black(), theme_latexfonts())) do
#     fig = Figure(resolution = (400, 250), backgroundcolor = (:black, 0.0), textcolor=:white, axiscolor=:white, color=:white)
#     ax = Axis(fig[1, 1], xlabel = L"\nu", ylabel = "fit error", backgroundcolor=(:black,0.0))
#     lw = 1.8
#     err_lw = 1.1
#     ms = 3
#     scatterlines!(ax, p_space[2], residuals[2], markersize = ms, color=:white)
#     save("YaoKivelsonOrientable_collapse_residuals.pdf", fig)
# end

fig, ax = scatter(p_space[2], residuals[2])
vlines!(ax, [sp.optimal_ps[2]], color=:red)
fig
sp.minimum
function scaling(x, L, nu)
    # return x.*L^(1/nu)
    return x
end
using ColorSchemes
blue = colorant"#175EC7"
red = colorant"#C81E5B"
colormap = ColorScheme(range(red, blue))
colors = reverse([get(colormap, x) for x in range(0.0, 1.0, 7)])
with_theme(merge(theme_black(), theme_latexfonts())) do
J = [params[i][1] for i in eachindex(params)]
K = [params[i][2] for i in eachindex(params)]
Kc = 0.63
x = (K .- Kc)

fig = Figure(resolution = (700, 400), backgroundcolor = (:black, 0.0), textcolor=:white, axiscolor=:white, color=:white)
ax = Axis(fig[1, 1], #=xlabel = L"(K-K_c) L^{1/\nu}",=# xlabel=L"(K-K_c)", ylabel = L"I(A:B:C)", backgroundcolor=(:black,0.0))
lw = 1.8
err_lw = 1.1
ms = 8
nu = 1.0

p_12 = scatterlines!(ax, scaling(x, 12, nu), I3_12, linewidth = lw, markersize = ms, label = "L = 12", color=colors[1])
errorbars!(ax, scaling(x, 12, nu), I3_12, I3_12_err, color=colors[1], linewidth=err_lw)
p_24 = scatterlines!(ax, scaling(x, 24, nu), I3_24, linewidth = lw, markersize = ms, label = "L = 24", color=colors[2])
errorbars!(ax, scaling(x, 24, nu), I3_24, I3_24_err, color=colors[2], linewidth=err_lw)
p_28 = scatterlines!(ax, scaling(x, 28, nu), I3_28, linewidth = lw, markersize = ms, label = "L = 28", color=colors[3])
errorbars!(ax, scaling(x, 28, nu), I3_28, I3_28_err, color=colors[3], linewidth=err_lw)
p_32 = scatterlines!(ax, scaling(x, 32, nu), I3_32, linewidth = lw, markersize = ms, label = "L = 32", color=colors[4])
errorbars!(ax, scaling(x, 32, nu), I3_32, I3_32_err, color=colors[4], linewidth=err_lw)
p_36 = scatterlines!(ax, scaling(x, 36, nu), I3_36, linewidth = lw, markersize = ms, label = "L = 36", color=colors[5])
errorbars!(ax, scaling(x, 36, nu), I3_36, I3_36_err, color=colors[5], linewidth=err_lw)
p_40 = scatterlines!(ax, scaling(x, 40, nu), I3_40, linewidth = lw, markersize = ms, label = "L = 40", color=colors[6])
errorbars!(ax, scaling(x, 40, nu), I3_40, I3_40_err, color=colors[6], linewidth=err_lw)
p_48 = scatterlines!(ax, scaling(x, 48, nu), I3_48, linewidth = lw, markersize = ms, label = "L = 48", color=colors[7])
errorbars!(ax, scaling(x, 48, nu), I3_48, I3_48_err, color=colors[7], linewidth=err_lw)
Legend(fig[1,2], reverse([p_12, p_24, p_28, p_32, p_36, p_40, p_48]), reverse(["L=12", "L=24", "L=28", "L=32", "L=36", "L=40", "L=48"]), bgcolor=(:black, 0.0))
fig
# save("YaoKivelsonOrientable_uncollapse.pdf", fig)
end