function sample_I3(circuit, thermalization_time, n_samples, sample_distance; filename = nothing)
    if filename != nothing
        if isfile(filename)
            return nothing
        end
    end
	ops = get_operators(circuit)
	I3 = 0.0
	state = initial_state(circuit)
	#thermalization
	for _ in 1:thermalization_time
		apply!(state, circuit, ops)
	end
	for s in 1:n_samples
		for _ in 1:sample_distance
			apply!(state, circuit, ops)
		end
		I3 += mean(tmi(state, circuit))
	end
	I3 = I3 / n_samples
	if filename != nothing
		jldsave("$(filename)"; I3 = I3)
	else
		return I3
	end
end
