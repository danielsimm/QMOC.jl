### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ efdf9154-56cc-11ee-172c-d97c3a1f54c1
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ cadc2aec-5e51-4da9-ad4e-4938fa359d36
# ╠═╡ show_logs = false
using LatticeCircuits

# ╔═╡ fb05c306-f352-4151-9fac-42553112a837
using PlutoUI

# ╔═╡ fc66e84b-1f09-40c0-8fbd-38546998a7cf
using CairoMakie

# ╔═╡ 74ef1f91-634f-48f9-9ab2-d1a5081b5a43
using CurveFit

# ╔═╡ 7c1b93f3-7f77-469d-9a2b-d78b33ffa57a
begin
	I3_1024 = evaluate(loadSimulation("PPcritical1024"), :I3)[2]
	I3_256 = evaluate(loadSimulation("PPcritical256"), :I3)[2]
	I3_512 = evaluate(loadSimulation("PPcritical512"), :I3)[2]
end;

# ╔═╡ 837a6cdf-348b-4984-a9b3-597a3790ae11
begin
	EE_1024 = evaluate(loadSimulation("PPcritical1024"), :SvN)[2]
	EE_256 = evaluate(loadSimulation("PPcritical256"), :SvN)[2]
	EE_512 = evaluate(loadSimulation("PPcritical512"), :SvN)[2]	
end;

# ╔═╡ d0f70fa1-fa60-4375-86c0-f16ff96bf4f1
xfunc(l) = (1/3)*log(sin(l*pi/64));

# ╔═╡ 85ecbdcf-0faf-4b15-a3f7-e3a89f2c4897
EE_x_axis = xfunc.(1:1:63);

# ╔═╡ 655ab63e-49c4-46d9-ba8e-17d2765dc5d7
params = loadSimulation("PPcritical1024").parameter_set;

# ╔═╡ 5c39da72-f370-4212-8d82-1ba9870a152b
@bind index PlutoUI.Slider(1:40)

# ╔═╡ 5549f662-6861-4615-9590-4eca1f0b81b6
begin
	fit_slopes1024 = zeros(40)
	fit_consts1024 = zeros(40)
	fit_slopes512 = zeros(40)
	fit_consts512 = zeros(40)
	ys1024 = [vec(EE_1024[index] .- EE_1024[index][32])[2:64] for index in 1:40]
	ys512 = [vec(EE_512[index])[2:64] .- EE_512[index][32] for index in 1:40]
	for index in eachindex(fit_slopes1024)
		fit1024 = linear_fit(EE_x_axis, ys1024[index])
		fit512 = linear_fit(EE_x_axis, ys512[index])
		fit_slopes1024[index] = fit1024[2]
		fit_consts1024[index] = fit1024[1]
		fit_slopes512[index] = fit512[2]
		fit_consts512[index] = fit512[1]
	end
end
		

# ╔═╡ 54cfb2bc-fb52-49c6-912c-6ad42ca5a1b7
begin
	fig = Figure()
	ax1 = Axis(fig[1,1], xlabel=L"l/L", ylabel=L"\Delta S_{vN}")
	ax2 = Axis(fig[2,1], xlabel=L"\frac{1}{3}\log(\sin(\pi \frac{l}{L}))", ylabel=L"\Delta S_{vN}")
	ax3 = Axis(fig[1,2], aspect = DataAspect())
	ax4 = Axis(fig[2,2], xlabel=L"p_z", ylabel=L"c_{eff}")
	y1024 = ys1024[index]
	y512 = ys512[index]
	pzs = [params[i][3] for i in 1:40]
	fitfunc1024(x) = fit_slopes1024[index]*x + fit_consts1024[index]
	fitfunc512(x) = fit_slopes512[index]*x + fit_consts512[index]
    lines!(ax3, [0.0, 1.0, 0.5, 0.0], [0.0, 0.0, sqrt(3)/2, 0.0], color = :black)
    param_x, param_y = LatticeCircuits.parametric_to_cartesian(params[index])
	scatter!(ax1, (1:1:63)./64, y1024, color=:red)
	scatter!(ax1, (1:1:63)./64, y512)
	scatter!(ax2, EE_x_axis, y1024, color=:red)
	scatter!(ax2, EE_x_axis, y512)
	lines!(ax2, -1:0, fitfunc1024.(-1:0), color=:red)
	lines!(ax2, -1:0, fitfunc512.(-1:0))
	lines!(ax4, pzs, fit_slopes1024, color=:red, label=L"L=1024")
	scatter!(ax4, pzs[index], fit_slopes1024[index], color=:red)
	lines!(ax4, pzs, fit_slopes512, label=L"L=512")
	hlines!(ax4, [(3*sqrt(3)/(pi))], color=:green, label=L"\frac{3 \sqrt{3}}{\pi}")
	axislegend(ax4, position=:lt)
	scatter!(ax4, pzs[index], fit_slopes512[index])
	scatter!(ax3, param_x, param_y, markersize = 15)
    hidedecorations!(ax3)
    hidespines!(ax3)
    text!(ax3, 0.5, sqrt(3)/2+0.05, text=L"p_z", fontsize=20, align=(:center, :center))
    text!(ax3, -0.05, -0.05, text=L"p_x", fontsize=20, align=(:center, :center))
    text!(ax3, 1.05, -0.05, text=L"p_y", fontsize=20, align=(:center, :center))
	fig
end

# ╔═╡ 0b0e26d3-588c-44e7-9445-2faf8238d226
SimulationArchive = LatticeCircuits.loadSimulationArchive()

# ╔═╡ 1d9ecf7b-df53-4418-a870-d67e45f21059
LatticeCircuits.loadSimulationArchive()

# ╔═╡ 9f6243f6-c0e0-4e14-a332-964c80be4122
evaluate(loadSimulation("PPPhasetransition128"), :I3)

# ╔═╡ Cell order:
# ╠═efdf9154-56cc-11ee-172c-d97c3a1f54c1
# ╠═cadc2aec-5e51-4da9-ad4e-4938fa359d36
# ╠═fb05c306-f352-4151-9fac-42553112a837
# ╠═7c1b93f3-7f77-469d-9a2b-d78b33ffa57a
# ╠═837a6cdf-348b-4984-a9b3-597a3790ae11
# ╠═d0f70fa1-fa60-4375-86c0-f16ff96bf4f1
# ╠═85ecbdcf-0faf-4b15-a3f7-e3a89f2c4897
# ╠═655ab63e-49c4-46d9-ba8e-17d2765dc5d7
# ╠═fc66e84b-1f09-40c0-8fbd-38546998a7cf
# ╠═74ef1f91-634f-48f9-9ab2-d1a5081b5a43
# ╠═5c39da72-f370-4212-8d82-1ba9870a152b
# ╟─5549f662-6861-4615-9590-4eca1f0b81b6
# ╟─54cfb2bc-fb52-49c6-912c-6ad42ca5a1b7
# ╠═0b0e26d3-588c-44e7-9445-2faf8238d226
# ╠═1d9ecf7b-df53-4418-a870-d67e45f21059
# ╠═9f6243f6-c0e0-4e14-a332-964c80be4122
