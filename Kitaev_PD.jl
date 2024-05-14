using QMOC
using Statistics
using CairoMakie
using GeometryBasics
using ColorSchemes
L = 32
res = 12
samples = 16
subsamples = 50

begin
    parameters = vec.(unique(QMOC.symmetry_data_extension(parameter_wedge(res))))
    circuits = [QMOC.KitaevCircuit(L, parameters[i]) for i in eachindex(parameters)]
    I3s = zeros(length(circuits))
    ops = QMOC.get_operators(circuits[1])
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
                I3 += mean(QMOC.tmi(state, c))
            end
            this_I3[s] = I3 / subsamples
        end
        I3s[i] = mean(this_I3)
    end
end


begin
    parameter_set = QMOC.parametric_to_cartesian.(parameters)
    unique_indices = unique(i -> parameter_set[i], eachindex(parameter_set))
    parameter_set = parameter_set[unique_indices]
    x = [parameter_set[i][1] for i in eachindex(parameter_set)]
    y = [parameter_set[i][2] for i in eachindex(parameter_set)]
    I3 = I3s[unique_indices]


    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = "non", titlesize = 25)


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

	fig
end


