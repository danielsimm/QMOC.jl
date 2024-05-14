using CairoMakie
using Colors
using ColorSchemes
using GeometryBasics
using QuantumClifford
using QMOC

const latticeblack = RGBf(28/256, 40/256, 51/256)
const latticered = RGBf(231/256, 76/256, 60/256)
const latticegreen = RGBf(33/256, 199/256, 135/256)
const latticeblue = RGBf(52/256, 152/256, 219/256)

"""
    parameter_legend(parameters)

    Plots all parameter sets in `parameters` on the parameter triangle, using `Makie.jl`.
"""
function parameter_legend(parameters)
    scheme = cgrad(:Dark2_7, rev=false)
    fig = Figure(resolution = (1000, 1000))
    ax = Axis(fig[1, 1], aspect = DataAspect())
    lines!(ax, [0.0, 1.0, 0.5, 0.0], [0.0, 0.0, sqrt(3)/2, 0.0], color = :black)
    x_values = zeros(length(parameters))
    y_values = zeros(length(parameters))
    z_values = zeros(length(parameters))
    colors = []
    for i in eachindex(z_values)
        x_values[i], y_values[i] = parametric_to_cartesian(parameters[i])
        z_values[i] = parameter_distance(parameters[i])
        push!(colors, scheme[z_values[i]])
    end
    # for i in eachindex(z_values)
    #     for j in eachindex(z_values)
    #         if i != j
    #             if z_values[i] == z_values[j]
    #                 lines!(ax, parameter_circle(z_values[i])[1], parameter_circle(z_values[i])[2], color = colors[i])
    #             end
    #         end
    #     end
    # end
    scatter!(ax, x_values, y_values, color = colors, markersize = 25)
    hidedecorations!(ax)
    hidespines!(ax)
    text!(ax, 0.5, sqrt(3)/2+0.05, text=L"p_z", fontsize=50, align=(:center, :center))
    text!(ax, -0.05, -0.05, text=L"p_x", fontsize=50, align=(:center, :center))
    text!(ax, 1.05, -0.05, text=L"p_y", fontsize=50, align=(:center, :center))
    return fig
end

function parameter_circles2(r)
    radius = r
    x_values1 = zeros(361)
    y_values1 = zeros(361)
    for i in 1:361
        x_values1[i] = radius*cos(i*pi/180) 
        y_values1[i] = radius*sin(i*pi/180)
    end 
    x_values2 = zeros(361)
    y_values2 = zeros(361)
    for i in 1:361
        x_values2[i] = radius*cos(i*pi/180) .+ 1.0
        y_values2[i] = radius*sin(i*pi/180)
    end
    x_values3 = zeros(361)
    y_values3 = zeros(361)
    for i in 1:361
        x_values3[i] = radius*cos(i*pi/180) .+ 0.5
        y_values3[i] = radius*sin(i*pi/180) .+ sqrt(3)/2
    end
    return x_values1, y_values1, x_values2, y_values2, x_values3, y_values3
end


### Honeycomb ###

function operatorqubits(operator::PauliOperator)
    qubits = []
    strings = []
    for i in eachindex(operator)
        if operator[i][1] || operator[i][2]
            push!(qubits, i)
            if operator[i][1]
                if operator[i][2]
                    push!(strings, "Y")
                else
                    push!(strings, "X")
                end
            else
                push!(strings, "Z")
            end
        end
    end
    return qubits, strings
end

function plot(lattice::HoneycombLattice, highlight = nothing; label=false)

    L = lattice.L
    sites = lattice.sites

    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1], aspect = DataAspect())

    for (i, site) in enumerate(sites)
        x, y = site.cartesianx, site.cartesiany
        xneighbour_x, xneighbour_y = sites[site.xneighbour].cartesianx, sites[site.xneighbour].cartesiany
        yneighbour_x, yneighbour_y = sites[site.yneighbour].cartesianx, sites[site.yneighbour].cartesiany
        zneighbour_x, zneighbour_y = sites[site.zneighbour].cartesianx, sites[site.zneighbour].cartesiany
        x_neighbour_distance = sqrt((x - xneighbour_x)^2 + (y - xneighbour_y)^2)
        y_neighbour_distance = sqrt((x - yneighbour_x)^2 + (y - yneighbour_y)^2)
        z_neighbour_distance = sqrt((x - zneighbour_x)^2 + (y - zneighbour_y)^2)

        if x_neighbour_distance < 1.5
            lines!(ax, [x, xneighbour_x], [y, xneighbour_y], color=latticered, linewidth=5)
        else
            lines!(ax, [x, xneighbour_x], [y, xneighbour_y], color=latticered, linewidth=1)
        end

        if y_neighbour_distance < 1.5
            lines!(ax, [x, yneighbour_x], [y, yneighbour_y], color=latticegreen, linewidth=5)
        else
           line = lines!(ax, [x, yneighbour_x], [y, yneighbour_y], color=latticegreen, linewidth=2, linestyle=:dash)
            translate!(line, 0, 0, -1)
        end

        if z_neighbour_distance < 1.5
            lines!(ax, [x, zneighbour_x], [y, zneighbour_y], color=latticeblue, linewidth=5)
        else
            line = lines!(ax,  [x, zneighbour_x], [y, zneighbour_y], color=latticeblue, linewidth=2, linestyle=:dash)
            translate!(line, 0, 0, -1)
        end

        point = scatter!(ax, x, y, markersize = 150/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(site.index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    if highlight !== nothing
        if highlight isa Int
            x = sites[highlight].cartesianx
            y = sites[highlight].cartesiany
            scatter!(ax, x, y, markersize = 200/L, color=(:yellow, 0.5))
        else
            for i in eachindex(highlight)
                x = sites[highlight[i]].cartesianx
                y = sites[highlight[i]].cartesiany
                scatter!(ax, x, y, markersize = 250/L, color=(:yellow, 0.5))
            end
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

plot(HoneycombLattice(8); label=true)


_HC_zrows(8)
_HC_yrows(8)
_HC_xrows(8)



_HC_all_partitions(8)

for partition in _HC_all_partitions(8)
    for i in eachindex(partition)
        fig = plot(HoneycombLattice(8), partition[i]; label=true)
        display(fig)
        sleep(1)
    end
end


function HC_plot_lattice_kekule(lattice::Vector{QMOC.HoneycombLatticeSite}, highlight=nothing; label=false)
    posterblue = colorant"#3498DB";
    postergreen = colorant"#46AD77";
    posterred = colorant"#E74C3C";
    posterblack = colorant"#1C2833";
    posterdarkblue = colorant"#2980B9";
    posterdark = colorant"#2C3E50";
    background = colorant"#ECF0F1";
    highlightcolor = colorant"#F29325";
    L = Int(sqrt(length(lattice) / 2))
    fig = Figure(resolution = (800, 800), backgroundcolor = (background, 0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = (background, 0.0))
    for i in 2:length(lattice)-1
        x = lattice[i].cartesianx
        y = lattice[i].cartesiany
        redneighbourx = lattice[lattice[i].redneighbour].cartesianx
        redneighboury = lattice[lattice[i].redneighbour].cartesiany
        greenneighbourx = lattice[lattice[i].greenneighbour].cartesianx
        greenneighboury = lattice[lattice[i].greenneighbour].cartesiany
        blueneighbourx = lattice[lattice[i].blueneighbour].cartesianx
        blueneighboury = lattice[lattice[i].blueneighbour].cartesiany
          
        if !(abs(x - redneighbourx) > 1.1) && lattice[i].redneighbour in 2:length(lattice)-1
            lines!(ax, [x, redneighbourx], [y, redneighboury], color=posterred, linewidth=4)
        end
        if !(abs(x - greenneighbourx) > 1.1) && lattice[i].greenneighbour in 2:length(lattice)-1
            lines!(ax, [x, greenneighbourx], [y, greenneighboury], color=postergreen, linewidth=4)
        end
        if !(abs(x - blueneighbourx) > 1.8) && lattice[i].blueneighbour in 2:length(lattice)-1
            lines!(ax, [x, blueneighbourx], [y, blueneighboury], color=posterblue, linewidth=4)
        end
        
        point = scatter!(ax, x, y, markersize = 150/L, color=posterblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(lattice[i].lindex)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    if highlight !== nothing
        for i in eachindex(highlight)
            x = lattice[highlight[i]].cartesianx
            y = lattice[highlight[i]].cartesiany
            scatter!(ax, x, y, markersize = 200/L, color=(:highlightcolor, 0.8))
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function DHC_plot_lattice(lattice::Vector{QMOC.DecoratedHoneycombLatticeSite}, highlight = nothing; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = background)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = background)
    L = Int(sqrt(div(length(lattice), 6)))
    for i in eachindex(lattice)
        x = lattice[i].cartesianx
        y = lattice[i].cartesiany
        xneighbourx = lattice[lattice[i].xneighbour].cartesianx
        xneighboury = lattice[lattice[i].xneighbour].cartesiany
        yneighbourx = lattice[lattice[i].yneighbour].cartesianx
        yneighboury = lattice[lattice[i].yneighbour].cartesiany
        zneighbourx = lattice[lattice[i].zneighbour].cartesianx
        zneighboury = lattice[lattice[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticered, linewidth=4)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticegreen, linewidth=4)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblue, linewidth=4)
        end
        point = scatter!(ax, x, y, markersize = 80/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(lattice[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    if highlight !== nothing
        for i in eachindex(highlight)
            x = lattice[highlight[i]].cartesianx
            y = lattice[highlight[i]].cartesiany
            scatter!(ax, x, y, markersize = 120/L, color=(:yellow, 1.0))
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function DHC_plot_lattice(lattice::Vector{QMOC.DecoratedHoneycombLatticeSite}, highlight::QuantumClifford.PauliOperator; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = background)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = background)
    L = Int(sqrt(div(length(lattice), 6)))
    for i in eachindex(lattice)
        x = lattice[i].cartesianx
        y = lattice[i].cartesiany
        xneighbourx = lattice[lattice[i].xneighbour].cartesianx
        xneighboury = lattice[lattice[i].xneighbour].cartesiany
        yneighbourx = lattice[lattice[i].yneighbour].cartesianx
        yneighboury = lattice[lattice[i].yneighbour].cartesiany
        zneighbourx = lattice[lattice[i].zneighbour].cartesianx
        zneighboury = lattice[lattice[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticered, linewidth=4.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticegreen, linewidth=4.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblue, linewidth=4.0)
        end
        point = scatter!(ax, x, y, markersize = 80/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(lattice[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    for i in eachindex(operatorqubits(highlight)[1])
        qubit = operatorqubits(highlight)[1][i]
        string = operatorqubits(highlight)[2][i]
        x = lattice[qubit].cartesianx
        y = lattice[qubit].cartesiany
        point = scatter!(ax, x, y, markersize = 120/L, color=(:yellow, 1.0))
        translate!(point,0,0,2)
        text = text!(x, y, text = string, color=:black, align=(:center, :center))
        translate!(text,0,0,3)
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function DHC_plot_lattice(lattice::Vector{QMOC.DecoratedHoneycombLatticeSite}, highlights::Vector{QuantumClifford.PauliOperator}; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = background)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = background)
    L = Int(sqrt(div(length(lattice), 6)))
    for i in eachindex(lattice)
        x = lattice[i].cartesianx
        y = lattice[i].cartesiany
        xneighbourx = lattice[lattice[i].xneighbour].cartesianx
        xneighboury = lattice[lattice[i].xneighbour].cartesiany
        yneighbourx = lattice[lattice[i].yneighbour].cartesianx
        yneighboury = lattice[lattice[i].yneighbour].cartesiany
        zneighbourx = lattice[lattice[i].zneighbour].cartesianx
        zneighboury = lattice[lattice[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticered, linewidth=4.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticegreen, linewidth=4.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblue, linewidth=4.0)
        end
        point = scatter!(ax, x, y, markersize = 80/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(lattice[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    for highlight in highlights
        for i in eachindex(operatorqubits(highlight)[1])
            qubit = operatorqubits(highlight)[1][i]
            string = operatorqubits(highlight)[2][i]
            x = lattice[qubit].cartesianx
            y = lattice[qubit].cartesiany
            point = scatter!(ax, x, y, markersize = 120/L, color=(:yellow, 1.0))
            translate!(point,0,0,2)
            text = text!(x, y, text = string, color=:black, align=(:center, :center))
            translate!(text,0,0,3)
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function plot_lattice(lattice::QMOC.DecoratedHoneycombLattice; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = (:white, 0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = (:white, 0.0))
    L = lattice.L
    sites = lattice.sites
    # draw sites
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        point = scatter!(ax, x, y, markersize = 80/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(sites[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    # draw lattice
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        xneighbourx = sites[sites[i].xneighbour].cartesianx
        xneighboury = sites[sites[i].xneighbour].cartesiany
        yneighbourx = sites[sites[i].yneighbour].cartesianx
        yneighboury = sites[sites[i].yneighbour].cartesiany
        zneighbourx = sites[sites[i].zneighbour].cartesianx
        zneighboury = sites[sites[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticeblack, linewidth=1.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticeblack, linewidth=1.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblack, linewidth=1.0)
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function plot_operators(lattice::QMOC.DecoratedHoneycombLattice, ops; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = (:white, 0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = (:white, 0.0))
    L = lattice.L
    sites = lattice.sites
    # draw sites
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        point = scatter!(ax, x, y, markersize = 80/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(sites[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    # draw lattice
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        xneighbourx = sites[sites[i].xneighbour].cartesianx
        xneighboury = sites[sites[i].xneighbour].cartesiany
        yneighbourx = sites[sites[i].yneighbour].cartesianx
        yneighboury = sites[sites[i].yneighbour].cartesiany
        zneighbourx = sites[sites[i].zneighbour].cartesianx
        zneighboury = sites[sites[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticeblack, linewidth=1.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticeblack, linewidth=1.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblack, linewidth=1.0)
        end
    end
    # draw operators
    for op in ops
        qubits = operatorqubits(op)[1]
        strings = operatorqubits(op)[2]
        if strings[1] == "X"
            color = latticered
        elseif strings[1] == "Y"
            color = latticegreen
        elseif strings[1] == "Z"
            color = latticeblue
        end
        x1 = sites[qubits[1]].cartesianx
        y1 = sites[qubits[1]].cartesiany
        x2 = sites[qubits[2]].cartesianx
        y2 = sites[qubits[2]].cartesiany
        if !(abs(x1 - x2) > 1.5)
            if !(abs(y1 - y2) > 1.5)
                lines!(ax, [x1, x2], [y1, y2], color=color, linewidth=4.0)
            end
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function plot_operators(lattice::QMOC.HoneycombLattice, ops; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(size = (800, 800), backgroundcolor = (:white, 0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = (:white, 0.0))
    L = lattice.L
    sites = lattice.sites
    # draw sites
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        point = scatter!(ax, x, y, markersize = 120/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(sites[i].index)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    # draw lattice
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        xneighbourx = sites[sites[i].xneighbour].cartesianx
        xneighboury = sites[sites[i].xneighbour].cartesiany
        yneighbourx = sites[sites[i].yneighbour].cartesianx
        yneighboury = sites[sites[i].yneighbour].cartesiany
        zneighbourx = sites[sites[i].zneighbour].cartesianx
        zneighboury = sites[sites[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticeblack, linewidth=2.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticeblack, linewidth = 2.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblack, linewidth = 2.0)
        end
    end
    # draw operators
    for op in ops
        qubits = operatorqubits(op)[1]
        strings = operatorqubits(op)[2]
        if strings[1] == "X" && strings[2] == "X"
            color = latticered
        elseif strings[1] == "Y" && strings[2] == "Y"
            color = latticegreen
        elseif strings[1] == "Z" && strings[2] == "Z"
            color = latticeblue
        else
            color = :orange
        end
        xs = [sites[qubits[i]].cartesianx for i in eachindex(qubits)]
        ys = [sites[qubits[i]].cartesiany for i in eachindex(qubits)]
        lines!(ax, xs, ys, color=color, linewidth=8.0)
        lines!(ax, [xs[1], xs[end]], [ys[1], ys[end]], color=color, linewidth=8.0)
        for i in eachindex(qubits)
            point = scatter!(ax, xs[i], ys[i], markersize = 120/L, color=latticeblack)
            translate!(point,0,0,3)
            text = text!(ax, xs[i], ys[i]; text = strings[i], color=:white, align=(:center, :center))
           translate!(text,0,0,4)
        end
        
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end
inds = [36, 38, 40, 59, 63, 55] .+18
plot_operators(QMOC.HoneycombLattice(4), [QMOC._HC_WilsonLoops(4)[1]]; label=true)
[QMOC._HC_WilsonLoops(4)...]
for i in eachindex(QMOC._HC_WilsonPlaquette_operators(4))
    fig = plot_operators(QMOC.HoneycombLattice(4), [QMOC._HC_WilsonPlaquette_operators(4)[i]]; label=true)
    display(fig)
    sleep(1.5)
end

all_ops = [QMOC._HC_XX_operators(4)..., QMOC._HC_YY_operators(4)..., QMOC._HC_ZZ_operators(4)...]
loop = QMOC._HC_WilsonLoops(4)[2]
ops = QMOC._HC_YY_operators(4)
for op in all_ops
    for loop in QMOC._HC_WilsonPlaquette_operators(4)
    if !(QuantumClifford.comm(op, loop) == 0)
        println("!!!")
    end
end
end

fig = plot_operators(QMOC.HoneycombLattice(4), [QMOC._HC_WilsonPlaquette_operators(4)[1]]; label=true)

QMOC._HC_WilsonPlaquette_operators(4)[1]


length(QMOC._HC_ZZ_operators(4))

function plot_operators_kekule(lattice::QMOC.HoneycombLattice, ops::Matrix{QuantumClifford.PauliOperator}; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    fig = Figure(resolution = (800, 800), backgroundcolor = (:white, 0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = (:white, 0.0))
    L = lattice.L
    sites = lattice.sites
    # draw sites
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        point = scatter!(ax, x, y, markersize = 120/L, color=latticeblack)
        translate!(point,0,0,1)
        if label
            text = text!(x, y, text = "$(sites[i].lindex)", color=:white, align=(:center, :center))
            translate!(text,0,0,2)
        end
    end
    # draw lattice
    for i in eachindex(sites)
        x = sites[i].cartesianx
        y = sites[i].cartesiany
        xneighbourx = sites[sites[i].xneighbour].cartesianx
        xneighboury = sites[sites[i].xneighbour].cartesiany
        yneighbourx = sites[sites[i].yneighbour].cartesianx
        yneighboury = sites[sites[i].yneighbour].cartesiany
        zneighbourx = sites[sites[i].zneighbour].cartesianx
        zneighboury = sites[sites[i].zneighbour].cartesiany
        lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticeblack, linewidth=2.0)
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticeblack, linewidth = 2.0)
        end
        if !(abs(y - zneighboury) > 1.5)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblack, linewidth = 2.0)
        end
    end
    # draw operators
    # for op in ops
    #     qubits = operatorqubits(op)[1]
    #     strings = operatorqubits(op)[2]
    #     if strings[1] == "X"
    #         color = latticered
    #     elseif strings[1] == "Y"
    #         color = latticegreen
    #     elseif strings[1] == "Z"
    #         color = latticeblue
    #     end
    #     x1 = sites[qubits[1]].cartesianx
    #     y1 = sites[qubits[1]].cartesiany
    #     x2 = sites[qubits[2]].cartesianx
    #     y2 = sites[qubits[2]].cartesiany
    #     if !(abs(x1 - x2) > 1.5)
    #         if !(abs(y1 - y2) > 1.5)
    #             lines!(ax, [x1, x2], [y1, y2], color=color, linewidth=8.0)
    #         end
    #     end
    # end
    for i in 1:L^2 # red connections
        op = ops[1, i]
        qubits = operatorqubits(op)[1]
        x1 = sites[qubits[1]].cartesianx
        y1 = sites[qubits[1]].cartesiany
        x2 = sites[qubits[2]].cartesianx
        y2 = sites[qubits[2]].cartesiany
        if !(abs(x1 - x2) > 1.5)
            if !(abs(y1 - y2) > 1.5)
                lines!(ax, [x1, x2], [y1, y2], color=latticered, linewidth=8.0)
            end
        end
    end
    for i in 1:L^2 # green connections
        op = ops[2, i]
        qubits = operatorqubits(op)[1]
        x1 = sites[qubits[1]].cartesianx
        y1 = sites[qubits[1]].cartesiany
        x2 = sites[qubits[2]].cartesianx
        y2 = sites[qubits[2]].cartesiany
        if !(abs(x1 - x2) > 1.5)
            if !(abs(y1 - y2) > 1.5)
                lines!(ax, [x1, x2], [y1, y2], color=latticegreen, linewidth=8.0)
            end
        end
    end
    for i in 1:L^2 # blue connections
        op = ops[3, i]
        qubits = operatorqubits(op)[1]
        x1 = sites[qubits[1]].cartesianx
        y1 = sites[qubits[1]].cartesiany
        x2 = sites[qubits[2]].cartesianx
        y2 = sites[qubits[2]].cartesiany
        if !(abs(x1 - x2) > 1.5)
            if !(abs(y1 - y2) > 1.5)
                lines!(ax, [x1, x2], [y1, y2], color=latticeblue, linewidth=8.0)
            end
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end


