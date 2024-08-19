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

function sample_full(circuit, thermalization_time, n_samples, sample_distance; filename = nothing)
	if filename != nothing
		if isfile(filename)
			return nothing
		end
	end
	ops = get_operators(circuit)
	I3 = 0.0
	EE = zeros(length(subsystem_labels(circuit)))
	state = initial_state(circuit)
	#thermalization
	for _ in 1:thermalization_time
		apply!(state, circuit, ops)
	end
	for s in 1:n_samples
		for _ in 1:sample_distance
			apply!(state, circuit, ops)
		end
		I3 += 0.0 #mean(tmi(state, circuit))
		EE .+= entropy(state, circuit)
	end
	I3 = I3 / n_samples
	EE = EE ./ n_samples
	if filename != nothing
		jldsave("$(filename)"; I3 = I3, EE = EE)
	else
		return I3, EE
	end
end

function sample_free_purification(circuit, timesteps; filename = nothing)
	if filename != nothing
		if isfile(filename)
			return nothing
		end
	end
	ops = get_operators(circuit)
	entropy = zeros(timesteps + 1)
	state = initial_mixed_state(circuit)
	entropy[1] = circuit.nqubits - QuantumClifford.rank(state)
	#thermalization
	for i in 1:timesteps
		apply!(state, circuit, ops)
		entropy[i+1] = circuit.nqubits - QuantumClifford.rank(state)
	end
	if filename != nothing
		jldsave("$(filename)"; entropy = entropy)
	else
		return entropy
	end
end

function sample_purification_cleanup(folder, timesteps)
	files = readdir(folder)
	to_read = filter(contains("idx"), files)
	idx_array = sort(unique(parse.(Int64, [split(split(string, "idx")[2], "_")[1] for string in to_read])))

	Threads.@threads for idx in idx_array
		this_files = filter(contains("idx$(idx)_"), to_read)
		entropies = zeros(length(this_files), timesteps + 1)
		for i in eachindex(this_files)
			file = this_files[i]
			entropies[i, :] = jldopen("$(folder)/$(file)")["entropy"]
			rm("$(folder)/$(file)")
		end
		means = mean(entropies, dims = 1)
		stds = std(entropies, dims = 1)
		errs = stds ./ sqrt(size(entropies, 1))
		writedlm("$(folder)/idx$(idx)_entropy_mean.txt", [means stds errs])
	end
end

function sample_full_purification(circuit, timesteps; filename = nothing)
	if filename != nothing
		if isfile(filename)
			return nothing
		end
	end
	ops = get_operators(circuit)
	entropy = zeros(timesteps + 1)
	state = QuantumClifford.MixedDestabilizer(zero(QuantumClifford.Stabilizer, 1, circuit.nqubits))
	entropy[1] = circuit.nqubits
	#thermalization
	for i in 1:timesteps
		apply_single!(state, circuit, ops)
		entropy[i+1] = circuit.nqubits - QuantumClifford.rank(state)
	end
	if filename != nothing
		jldsave("$(filename)"; entropy = entropy)
	else
		return entropy
	end
end