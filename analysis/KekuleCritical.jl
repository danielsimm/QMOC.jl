using QMOC
using CairoMakie
using ScalingCollapse
using ColorSchemes

begin
	params, I3_12 = evaluate("KekuleCritical12", :I3)
	I3_24 = evaluate("KekuleCritical24", :I3)[2]
	I3_36 = evaluate("KekuleCritical36", :I3)[2]
	I3_48 = evaluate("KekuleCritical48", :I3)[2]
	I3_60 = evaluate("KekuleCritical60", :I3)[2]
	I3_72 = evaluate("KekuleCritical72", :I3)[2]
	I3_84 = evaluate("KekuleCritical84", :I3)[2]
	I3_12_err = evaluate("KekuleCritical12", :ΔI3)[2]
	I3_24_err = evaluate("KekuleCritical24", :ΔI3)[2]
	I3_36_err = evaluate("KekuleCritical36", :ΔI3)[2]
	I3_48_err = evaluate("KekuleCritical48", :ΔI3)[2]
	I3_60_err = evaluate("KekuleCritical60", :ΔI3)[2]
	I3_72_err = evaluate("KekuleCritical72", :ΔI3)[2]
	I3_84_err = evaluate("KekuleCritical84", :ΔI3)[2]

	I3_12 = Float64.(I3_12)
	I3_24 = Float64.(I3_24)
	I3_36 = Float64.(I3_36)
	I3_48 = Float64.(I3_48)
	I3_60 = Float64.(I3_60)
	I3_72 = Float64.(I3_72)
	I3_84 = Float64.(I3_84)
	I3_12_err = Float64.(I3_12_err) ./ sqrt(200)
	I3_24_err = Float64.(I3_24_err) ./ sqrt(200)
	I3_36_err = Float64.(I3_36_err) ./ sqrt(200)
	I3_48_err = Float64.(I3_48_err) ./ sqrt(200)
	I3_60_err = Float64.(I3_60_err) ./ sqrt(100)
	I3_72_err = Float64.(I3_72_err) ./ sqrt(100)
	I3_84_err = Float64.(I3_84_err) ./ sqrt(100)
end


begin
	x = Float64.([params[i][1] for i in eachindex(params)])
	data = [I3_12 I3_24 I3_36 I3_48 I3_60 I3_72 I3_84]
	err = [I3_12_err I3_24_err I3_36_err I3_48_err I3_60_err I3_72_err I3_84_err]
	perm = sortperm(x)
	x = x[perm]
	data = data[perm, :]
	err = err[perm, :]
end

data_12 = [x data[:, 1] err[:, 1]]
data_24 = [x data[:, 2] err[:, 2]]
data_36 = [x data[:, 3] err[:, 3]]
data_48 = [x data[:, 4] err[:, 4]]
data_60 = [x data[:, 5] err[:, 5]]
data_72 = [x data[:, 6] err[:, 6]]
data_84 = [x data[:, 7] err[:, 7]]

cd("../Master/plots/kekule/critical/data")
pwd()
using DelimitedFiles
writedlm("KekuleCritical12.txt", data_12)
writedlm("KekuleCritical24.txt", data_24)
writedlm("KekuleCritical36.txt", data_36)
writedlm("KekuleCritical48.txt", data_48)
writedlm("KekuleCritical60.txt", data_60)
writedlm("KekuleCritical72.txt", data_72)
writedlm("KekuleCritical84.txt", data_84)


sp = ScalingProblem(x, data, err, [12, 24, 36, 48, 60, 72, 84];
	sf = ScalingFunction(:x),
	dx = [0, Inf],
	p_space = [0.6:0.01:0.7, 0.9:0.01:1.1],
	#verbose=true,
	error = true,
)


function scaling(x, L, nu)
	return x .* L^(1 / nu)
end
using ColorSchemes
blue = colorant"#175EC7"
red = colorant"#C81E5B"
colormap = ColorScheme(range(red, blue))
colors = reverse([get(colormap, x) for x in range(0.0, 1.0, 7)])
with_theme(merge(theme_black(), theme_latexfonts())) do
	x = Float64.([params[i][1] for i in eachindex(params)])

	x = x .- 0.6247

	fig = Figure(resolution = (700, 400), backgroundcolor = (:black, 1.0), textcolor = :white, axiscolor = :white, color = :white)
	ax = Axis(fig[1, 1], xlabel = L"(p-p_c)L^{1/\nu}", ylabel = L"I(A:B:C)", backgroundcolor = (:black, 1.0))
	lw = 1.8
	err_lw = 1.1
	ms = 8
	nu = 1.047

	p_12 = scatter!(ax, scaling(x, 12, nu), I3_12, linewidth = lw, markersize = ms, label = "L = 12", color = (colors[1]))
	errorbars!(ax, scaling(x, 12, nu), I3_12, I3_12_err, color = colors[1], linewidth = err_lw)
	p_24 = scatter!(ax, scaling(x, 24, nu), I3_24, linewidth = lw, markersize = ms, label = "L = 24", color = (colors[2]))
	errorbars!(ax, scaling(x, 24, nu), I3_24, I3_24_err, color = colors[2], linewidth = err_lw)
	p_36 = scatter!(ax, scaling(x, 36, nu), I3_36, linewidth = lw, markersize = ms, label = "L = 36", color = (colors[3]))
	errorbars!(ax, scaling(x, 36, nu), I3_36, I3_36_err, color = colors[3], linewidth = err_lw)
	p_48 = scatter!(ax, scaling(x, 48, nu), I3_48, linewidth = lw, markersize = ms, label = "L = 48", color = (colors[4]))
	errorbars!(ax, scaling(x, 48, nu), I3_48, I3_48_err, color = colors[4], linewidth = err_lw)
	p_60 = scatter!(ax, scaling(x, 60, nu), I3_60, linewidth = lw, markersize = ms, label = "L = 60", color = (colors[5]))
	errorbars!(ax, scaling(x, 60, nu), I3_60, I3_60_err, color = colors[5], linewidth = err_lw)
	p_72 = scatter!(ax, scaling(x, 72, nu), I3_72, linewidth = lw, markersize = ms, label = "L = 72", color = (colors[6]))
	errorbars!(ax, scaling(x, 72, nu), I3_72, I3_72_err, color = colors[6], linewidth = err_lw)
	p_84 = scatter!(ax, scaling(x, 84, nu), I3_84, linewidth = lw, markersize = ms, label = "L = 84", color = (colors[7]))
	errorbars!(ax, scaling(x, 84, nu), I3_84, I3_84_err, color = colors[7], linewidth = err_lw)
	Legend(fig[1, 2], reverse([p_12, p_24, p_36, p_48, p_60, p_72, p_84]), reverse(["L=12", "L=24", "L=36", "L=48", "L=60", "L=72", "L=84"]), bgcolor = (:black, 0.0))

	vlines!(ax, [-0.5, 2.0], color = :white, linestyle = :dash, linewidth = 0.5)
	xlims!(ax, -5.0, 5.0)
	fig
	#save("Kekule_scaling.pdf", fig)
end
