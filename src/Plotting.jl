using CairoMakie
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

function HC_plot_lattice(lattice::Vector{HCSite}, highlight = nothing; label=false)
    background = RGBf(236/256, 240/256, 241/256)
    latticeblack = RGBf(28/256, 40/256, 51/256)
    latticered = RGBf(231/256, 76/256, 60/256)
    latticegreen = RGBf(33/256, 199/256, 135/256)
    latticeblue = RGBf(52/256, 152/256, 219/256)
    L = Int(sqrt(length(lattice) / 2))
    fig = Figure(resolution = (800, 800), backgroundcolor = background)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = background)
    for i in eachindex(lattice)
        x = lattice[i].cartesianx
        y = lattice[i].cartesiany
        xneighbourx = lattice[lattice[i].xneighbour].cartesianx
        xneighboury = lattice[lattice[i].xneighbour].cartesiany
        yneighbourx = lattice[lattice[i].yneighbour].cartesianx
        yneighboury = lattice[lattice[i].yneighbour].cartesiany
        zneighbourx = lattice[lattice[i].zneighbour].cartesianx
        zneighboury = lattice[lattice[i].zneighbour].cartesiany
        if !(abs(x - xneighbourx) > 1.1)
            lines!(ax, [x, xneighbourx], [y, xneighboury], color=latticeblack, linewidth=4)
        end
        if !(abs(x - yneighbourx) > 1.1)
            lines!(ax, [x, yneighbourx], [y, yneighboury], color=latticeblack, linewidth=4)
        end
        if !(abs(y - zneighboury) > 1.8)
            lines!(ax, [x, zneighbourx], [y, zneighboury], color=latticeblack, linewidth=4)
        end
        point = scatter!(ax, x, y, markersize = 150/L, color=latticeblack)
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
            scatter!(ax, x, y, markersize = 200/L, color=(:yellow, 0.5))
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function HC_plot_lattice_kekule(lattice::Vector{HCSite}, highlight=nothing; label=false)
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


function DHC_plot_lattice(lattice, highlight = nothing; label=false)
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
            text = text!(x, y, text = "$(lattice[i].lindex)", color=:white, align=(:center, :center))
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

function DHC_plot_lattice(lattice, highlight::QuantumClifford.PauliOperator; label=false)
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
            text = text!(x, y, text = "$(lattice[i].lindex)", color=:white, align=(:center, :center))
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

function DHC_plot_lattice(lattice, highlights::Vector{QuantumClifford.PauliOperator}; label=false)
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
            text = text!(x, y, text = "$(lattice[i].lindex)", color=:white, align=(:center, :center))
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