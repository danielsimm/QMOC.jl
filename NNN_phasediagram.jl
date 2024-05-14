using QMOC
using Statistics
using CairoMakie
using GeometryBasics
using ColorSchemes
L = 24
res = 12
samples = 16
subsamples = 100


Js = [i // 48 for i in 0:24]

function parameter_third(resolution)
    parameters = []
    starting_points_1 = parameter_line(:pz, [0, 1//2, 1//2], resolution)
    ending_points_1 = parameter_line(:pz, :center, resolution)
    starting_points_2 = parameter_line(:pz, [1//2, 0, 1//2], resolution)
    ending_points_2 = parameter_line(:pz, :center, resolution)
    for i in eachindex(starting_points_1)
        push!(parameters, parameter_line(starting_points_1[i], ending_points_1[i], resolution))
        push!(parameters, parameter_line(starting_points_2[i], ending_points_2[i], resolution))
    end
    return unique(reduce(vcat, parameters))
end

function Z3_symmetry_expand(params, observable)
    expanded_params = []
    expanded_observable = []
    for i in eachindex(params)
        p = params[i]
        o = observable[i]
        push!(expanded_params, p)
        push!(expanded_observable, o)
        push!(expanded_params, [p[2], p[3], p[1]])
        push!(expanded_observable, o)
        push!(expanded_params, [p[3], p[1], p[2]])
        push!(expanded_observable, o)
    end
    unique_indices = unique(i -> expanded_params[i], eachindex(expanded_params))
    return expanded_params[unique_indices], expanded_observable[unique_indices]
end

unique(QMOC.symmetry_data_extension(parameter_wedge(res)))





for (i, J) in enumerate(Js)


	wedge = parameter_wedge(res)
	parameters = [zeros(6) for _ in eachindex(wedge)]
	for i in eachindex(parameters)
		parameters[i][1:3] = wedge[i] .* (1 - J)
		parameters[i][4:6] .= (J) // 3
	end
	@assert all(p -> sum(p) ≈ 1, parameters)



	circuits = [QMOC.KitaevNNNCircuit(L, parameters[i]) for i in eachindex(parameters)]
	I3s = zeros(length(circuits))
	ops = QMOC.get_operators(circuits[1])
	QMOC.tmi(QMOC.initial_state(circuits[1]), circuits[1])
	length(QMOC.initial_state(circuits[1]))
	for (i, c) in enumerate(circuits)
		this_I3 = zeros(samples)
		Threads.@threads for s in 1:samples
			state = QMOC.initial_state(c)
			#thermalization
			for _ in 1:3*L
				QMOC.apply!(state, c, ops)
			end
			I3 = 0.0
			for _ in 1:subsamples
				QMOC.apply!(state, c, ops)
				I3 += QMOC.tmi(state, c)
			end
			this_I3[s] = I3 / subsamples
		end
		I3s[i] = mean(this_I3)
	end



	parameter_set = parameter_wedge(res)
	parameter_set, I3 = QMOC.symmetry_data_extension(parameter_set, I3s)
	parameter_set = QMOC.parametric_to_cartesian.(parameter_set)
	unique_indices = unique(i -> parameter_set[i], eachindex(parameter_set))
	parameter_set = parameter_set[unique_indices]
	x = [parameter_set[i][1] for i in eachindex(parameter_set)]
	y = [parameter_set[i][2] for i in eachindex(parameter_set)]
	I3 = I3[unique_indices]


	fig = Figure()
	ax = Axis(fig[1, 1], aspect = DataAspect(), title = L"$p_J= \frac{%$(i-1)}{40} $, $L=24$, 133 datapoints @ 1600 samples", titlesize = 25)


	tricontourf!(ax, x, y, Float64.(I3); show_generators = false, strokewidth = 0.0, colormap = reverse(ColorSchemes.tol_sunset), levels = -1:0.05:1, fxaa = true)
	hidedecorations!(ax)

	hidespines!(ax)

	xlims!(ax, -0.11, 1.11)
	ylims!(ax, -0.15, sqrt(3) / 2 + 0.09)
	triangle = Point2f[(0, 0), (0.5, sqrt(3) / 2), (1, 0), (0, 0)]
	poly!(ax, Polygon(Point2f[(-0.2, -0.2), (1.2, -0.2), (1.2, 1.2), (-0.2, 1.2)], [triangle]), color = (:white, 1.0))
	lines!(ax, triangle, color = :white, linewidth = 0)
	text!(ax, 0.5, sqrt(3) / 2 + 0.06, text = L"p_z", fontsize = 20, align = (:center, :center), color = :black)
	text!(ax, -0.05, -0.03, text = L"p_x", align = (:center, :center), color = :black, fontsize = 20)
	text!(ax, 1.05, -0.03, text = L"p_y", align = (:center, :center), color = :black, fontsize = 20)
	text!(ax, 0.5, -0.07, text = L"I_3", align = (:center, :center), color = :black, fontsize = 25)
	Colorbar(fig[1, 2], tellheight = false, colorrange = (-1, 1), colormap = reverse(ColorSchemes.tol_sunset))
	save("NNN_phasediagram_sweep_$(i).png", fig, px_per_unit = 2)
	fig
end


