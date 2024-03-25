function thermal_state(c::AbstractCircuit, time=3*c.size)
    stabilizer = initial_state(c)
    operators = get_operators(c)
    for t in 1:time
        apply!(stabilizer, c, operators)
    end
    return stabilizer
end