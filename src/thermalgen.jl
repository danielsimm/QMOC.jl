

function thermal_state(c::AbstractCircuit, time=3*c.size)
    stabilizer = initial_state(c)
    operators = get_operators(c)
    for t in 1:time
        apply!(stabilizer, c, operators)
    end
    return stabilizer
end

function write(state, c, name)
    if !isdir("$(hash(c))")
        mkdir("$(hash(c))")
    end
    jldsave("$(hash(c))/$(name).jld2", state=state)
end

function read(c, name)
    return jldopen("$(hash(c))/$(name).jld2") do file
        file["state"]
    end
end
