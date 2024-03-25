using QMOC
using Statistics
using CairoMakie
using GeometryBasics
if !isdir("temptestdata")
    mkdir("temptestdata")
end

cd("temptestdata")

# generate parameter set
parameter_set = parameter_wedge(15)
circuits = [QMOC.KitaevCircuit(8, params) for params in parameter_set]
n_states = 2
n_samples = 200

I3 = zeros(length(circuits))

Threads.@threads for i in eachindex(circuits)
    c = circuits[i]
    operators = QMOC.get_operators(c)
    this_I3 = zeros(n_states)
    for st in 1:n_states
        state = QMOC.thermal_state(c)
        for _ in 1:n_samples
            QMOC.apply!(state, c, operators)
            this_I3[st] += QMOC.tmi(state, c)
        end
        this_I3[st] /= n_samples
    end
    I3[i] = mean(this_I3)
end

# Threads.@threads for c in circuits
#     for i in 1:n_states
#         if !isfile("temptestdata/$(hash(c))_$(i).jld2")
#             stabilizer = QMOC.thermal_state(c)
#             QMOC.write(stabilizer, c, i)
#         end
#     end
# end


# Threads.@threads for i in eachindex(circuits)
#     c = circuits[i]
#     operators = QMOC.get_operators(c)
#     this_I3 = zeros(n_states)
#     for st in 1:n_states
#         state = QMOC.read(c, st)
#         for _ in 1:n_samples
#             QMOC.apply!(state, c, operators)
#             this_I3[st] += QMOC.tmi(state, c)
#         end
#         this_I3[st] /= n_samples
#     end
#     I3[i] = mean(this_I3)
# end

# generate plot
parameter_set, I3 = QMOC.symmetry_data_extension(parameter_set, I3)
parameter_set = QMOC.parametric_to_cartesian.(parameter_set)
unique_indices = unique(i -> parameter_set[i], eachindex(parameter_set))
parameter_set = parameter_set[unique_indices]
x = [parameter_set[i][1] for i in eachindex(parameter_set)]
y = [parameter_set[i][2] for i in eachindex(parameter_set)]
I3 = I3[unique_indices]

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())

scatter!(ax, x, y, color = I3.+1,markersize = 7, colormap = :viridis)
hidedecorations!(ax)

hidespines!(ax)

xlims!(ax, -0.11, 1.11)
ylims!(ax, -0.15, sqrt(3) / 2 + 0.09)
triangle = Point2f[(0, 0), (0.5, sqrt(3) / 2), (1, 0), (0, 0)]
poly!(ax, Polygon(Point2f[(-0.2, -0.2), (1.2, -0.2), (1.2, 1.2), (-0.2, 1.2)], [triangle]), color = (:white, 1.0))
lines!(ax, triangle, color = :black, linewidth = 1)
text!(ax, 0.5, sqrt(3) / 2 + 0.06, text = L"p_b", align = (:center, :center), color = :black)
text!(ax, -0.05, -0.03, text = L"p_r", align = (:center, :center), color = :black)
text!(ax, 1.05, -0.03, text = L"p_g", align = (:center, :center), color = :black)
text!(ax, 0.5, -0.07, text = L"I_3", align = (:center, :center), color = :black)
cd("..")
save("testplot.png", fig)
rm("temptestdata", recursive=true, force=true)
