using QMOC
using CairoMakie
using GeometryBasics
using ColorSchemes
set_theme!()
CairoMakie.activate!(type = "pdf", pt_per_unit = 1.0)


with_theme(theme_latexfonts()) do

blue = colorant"#175EC7"
red = colorant"#C81E5B"
params, I3 = evaluate("KekuleWedge20_L24", :I3)
params, I3 = QMOC.symmetry_data_extension(params, I3)
unique_indices = unique(i -> params[i], eachindex(params))
params = QMOC.parametric_to_cartesian.(params)
critical_px = 0.625
critical_pz = 0.096
transition_point1 = [critical_px, (1-critical_px)/2, (1-critical_px)/2]
transition_point2 = [(1-critical_pz)/2, (1-critical_pz)/2, critical_pz]
transition1_xy = Point2f(QMOC.parametric_to_cartesian(transition_point1))
transition2_xy = Point2f(QMOC.parametric_to_cartesian(transition_point2))
params = params[unique_indices]
I3 = I3[unique_indices]
I3 = Float64.(I3)
x = [params[i][1] for i in eachindex(params)]
y = [params[i][2] for i in eachindex(params)]
colormap = ColorScheme(range(red, blue))
fig = Figure(resolution = (1000, 1000), backgroundcolor = (:black, 0.0))
ax = Axis(fig[1,1],aspect=1)
xlims!(ax, -0.1, 1.1)
ylims!(ax, -0.1, sqrt(3)/2+0.1)
voronoiplot!(ax, x, y, I3, show_generators=false, strokewidth=0.01, colormap=colormap)
triangle = Point2f[(0,0),(0.5,sqrt(3)/2),(1,0),(0,0)]
poly!(Polygon(Point2f[(-0.2, -0.2), (1.2, -0.2), (1.2, 1.2), (-0.2, 1.2)], [triangle]), color = (:black, 1.0))
lines!(ax,triangle,color=:white, linewidth=1)
text!(ax, 0.5, sqrt(3)/2+0.05, text=L"p_b", fontsize=40, align=(:center, :center), color=:white)
text!(ax, -0.05, -0.05, text=L"p_r", fontsize=40, align=(:center, :center), color=:white)
text!(ax, 1.05, -0.05, text=L"p_g", fontsize=40, align=(:center, :center), color=:white)
Colorbar(fig[2,1], limits=(-1,1),colormap = colormap, labelsize=40, ticklabelsize=30, ticklabelcolor=:white, tickcolor=:white, vertical = false, width=700, height=20, alignmode=Outside(), tellwidth=false, labelpadding=-100)
Label(fig[3,1], L"I(A:B:C)", fontsize=40, color=:white, tellwidth=false)
scatter!(ax, [transition1_xy, transition2_xy], color=:orange, markersize=20)
hidespines!(ax)
hidedecorations!(ax)
fig
#save("KekulePhaseDiagram.pdf", fig)
end
params, I3 = evaluate("KekuleWedge20_L24", :I3)
params, I3 = QMOC.symmetry_data_extension(params, I3)
unique_indices = unique(i -> params[i], eachindex(params))
# params = QMOC.parametric_to_cartesian.(params)
critical_px = 0.625
critical_pz = 0.096
transition_point1 = [critical_px, (1-critical_px)/2, (1-critical_px)/2]
transition_point2 = [(1-critical_pz)/2, (1-critical_pz)/2, critical_pz]
transition1_xy = Point2f(QMOC.parametric_to_cartesian(transition_point1))
transition2_xy = Point2f(QMOC.parametric_to_cartesian(transition_point2))
params = params[unique_indices]
I3 = I3[unique_indices]
I3 = Float64.(I3)
x = [params[i][1] for i in eachindex(params)]
y = [params[i][2] for i in eachindex(params)]
params
data = [x y I3]

data = reduce(vcat, ([[params[i][1] params[i][2] params[i][3] I3[i]] for i in eachindex(params)]))

using DelimitedFiles
writedlm("../Master/plots/kekule/phase diagram/KekulePhaseDiagramData.dat", data)

append!(names, data)

# count unique x values
xvals = unique(x)
unique(y)

delaun