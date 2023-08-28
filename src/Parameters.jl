"""
   parameter_line(pointA, pointB, steps)

   Returns a list of rational parameter sets (px, py, pz) on a line between pointA and pointB with steps number of points.
"""
function parameter_line(pointA, pointB, steps)
    parameters = []
    if pointA == :center
        pointA = [1//3, 1//3, 1//3]
    elseif pointA == :px
        pointA = [1, 0, 0]
    elseif pointA == :py
        pointA = [0, 1, 0]
    elseif pointA == :pz
        pointA = [0, 0, 1]
    end
    if pointB == :center
        pointB = [1//3, 1//3, 1//3]
    elseif pointB == :px
        pointB = [1, 0, 0]
    elseif pointB == :py
        pointB = [0, 1, 0]
    elseif pointB == :pz
        pointB = [0, 0, 1]
    end
    step = -(pointA .- pointB)//(steps-1)
    for i in 1:steps
        px = pointA[1] + step[1]*(i-1)
        py = pointA[2] + step[2]*(i-1)
        pz = pointA[3] + step[3]*(i-1)
        push!(parameters, [px, py, pz])
    end
    return parameters
end

"""
    standard_lines(resolution)

    Returns a list of rational parameter sets (px, py, pz) on the two lines between the isotropic point and the px=1 point, px=py point.
"""
function standard_lines(resolution)
    parameters = []
    push!(parameters, parameter_line(:px, :center, resolution))
    push!(parameters, parameter_line([1//2, 1//2, 0], :center, resolution))
    return reduce(vcat, parameters)
end

"""
    parameter_wedge(resolution, mode=:radial)

    Returns a list of rational parameter sets (px, py, pz) filling the wedge between the isotropic point and the px=1 point, px=py point. Mode can be either :radial or :linear.
"""
function parameter_wedge(resolution, mode=:radial)
    starting_points = parameter_line(:px, [1//2, 1//2, 0], resolution)
    ending_points = parameter_line(:px, :center, resolution)
    parameters = []
    if mode == :radial
        for i in eachindex(starting_points)
            push!(parameters, parameter_line(starting_points[i], :center, resolution))
        end
    elseif mode == :linear
        for i in eachindex(starting_points)
            if starting_points[i] == ending_points[i]
                push!(parameters, [starting_points[i]])
            else
                push!(parameters, parameter_line(starting_points[i], ending_points[i], resolution))
            end
        end
    end
    return reduce(vcat, parameters)
end

"""
    parameter_distance(param)

    Returns the distance of a parameter set (px, py, pz) from the isotropic point in the parameter space.
"""
function parameter_distance(param)
    diff = param .- [1//3, 1//3, 1//3]
    return sqrt(sum(diff.^2))/sqrt(2//3)
end

"""
    parametric_to_cartesian(params)

    Returns the cartesian coordinates of a parameter set (px, py, pz).
"""
function parametric_to_cartesian(params)
    a = params[1]
    b = params[2]
    c = params[3]
    return [0.5*(2*b+c)/(a+b+c), (sqrt(3)*c)/(2*(a+b+c))]
end

"""
    parameter_circle(radius)

    Returns the cartesian coordinates of a circle with radius `radius` around the isometric point.
"""
function parameter_circle(radius)
    radius = radius/sqrt(3)
    x_values = zeros(361)
    y_values = zeros(361)
    for i in 1:361
        x_values[i] = radius*cos(i*pi/180) .+ 0.5
        y_values[i] = radius*sin(i*pi/180) .+ sqrt(3)/6
    end
    return x_values, y_values
end

# px critical: parameter_line([47//57, 5//57, 5//57], [29//57, 14//57, 14//57], 20)
