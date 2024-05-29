function starttime()
	println()
	println()
	println()
	println()
	println("---------------------------------------------")
	println("start:", now())
	println("---------------------------------------------")
	println()
	println()
	println()
	println()
end
function finaltime()
	println()
	println()
	println()
	println()
	println("---------------------------------------------")
	println("complete:", now())
	println("---------------------------------------------")
	println()
	println()
	println()
	println()
end

function distribute_work(circuits, n_trajectories)
	work = []
	for i in 1:n_trajectories
		idxs = randperm(length(circuits))
		for j in 1:length(circuits)
			push!(work, (i, idxs[j]))
		end
	end
	return work # (trajectory_id, circuit_id)
end

function meta_init(
	circuits::Vector{T} where T <: AbstractCircuit,
	n_trajectories::Int64,
	thermalization_time::Int64,
	n_samples::Int64,
	sample_distance::Int64,
	outputname::String,
)
	if !(isdir("cluster/$(outputname)"))
		mkdir("cluster/$(outputname)")
	end
	metadata = Dict(
		"circuits" => circuits,
		"n_trajectories" => n_trajectories,
		"thermalization_time" => thermalization_time,
		"n_samples" => n_samples,
		"sample_distance" => sample_distance,
		"outputname" => outputname,
	)
	jldsave("cluster/$(outputname)/metadata.jld2"; metadata = metadata)
end

function meta_init(meta::Dict, outputname::String)
	if !(isdir("cluster/$(outputname)"))
		mkdir("cluster/$(outputname)")
	end
	jldsave("cluster/$(outputname)/metadata.jld2"; metadata = meta)
end

function mpi_sample_I3(
	circuits::Vector{T} where T <: AbstractCircuit,
	n_trajectories::Int64,
	thermalization_time::Int64,
	n_samples::Int64,
	sample_distance::Int64,
	outputname::String,
)

	### MPI startup ###
	MPI.Init()
	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	nworkers = MPI.Comm_size(comm)
	root = 0
	MPI.Barrier(comm)
	###################

	### distribute work ###
	if rank == root
		meta_init(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, outputname)
		work = distribute_work(circuits, n_trajectories)
		part = [work[i:nworkers:end] for i in 1:nworkers]
		for i in 2:nworkers
			MPI.send(part[i], comm; dest = (i - 1))
		end
		todo = part[1]
		starttime()
		#println("rank $rank | has indices $(todo)")
	else
		todo = MPI.recv(comm) # recieve indices
		#println("rank $rank | has indices $(todo)")
	end
	MPI.Barrier(comm)
	########################

	if Threads.nthreads() > length(todo)
		println("rank $rank | Warning: more threads than work")
	end
	if Threads.nthreads() <= length(todo)
		println("rank $rank | running with $(Threads.nthreads()) threads | $(length(todo)) trajectories to sample...")
	end

	### do work ###
	Threads.@threads for this_work in todo
		trajectory_id, circuit_id = this_work
		circuit = circuits[circuit_id]
		filename = "cluster/$(outputname)/idx$(circuit_id)_$(trajectory_id).jld2"
		sample_I3(circuit, thermalization_time, n_samples, sample_distance; filename = filename)
		println("rank $rank | trajectory $(trajectory_id) of circuit $(circuit_id) -- done")
	end
	###############

	MPI.Barrier(comm)

	if rank == root
		folder = "cluster/$(outputname)"
		files = readdir(folder)
		to_read = filter(contains("idx"), files)
		idx_array = sort(unique(parse.(Int64, [split(split(string, "idx")[2], "_")[1] for string in to_read])))

		Threads.@threads for idx in idx_array
			this_files = filter(contains("idx$(idx)_"), to_read)
			out = zeros(length(this_files))
			for i in eachindex(this_files)
				file = this_files[i]
				out[i] = jldopen("$(folder)/$(file)")["I3"]
				rm("$(folder)/$(file)")
			end
			writedlm("$(folder)/idx$(idx).txt", out)
		end
		finaltime()
		MPI.Finalize()
	end
end

function mpi_sample_full(
	circuits::Vector{T} where T <: AbstractCircuit,
	n_trajectories::Int64,
	thermalization_time::Int64,
	n_samples::Int64,
	sample_distance::Int64,
	outputname::String,
)

	### MPI startup ###
	MPI.Init()
	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	nworkers = MPI.Comm_size(comm)
	root = 0
	MPI.Barrier(comm)
	###################

	### distribute work ###
	if rank == root
		meta_init(circuits, n_trajectories, thermalization_time, n_samples, sample_distance, outputname)
		work = distribute_work(circuits, n_trajectories)
		part = [work[i:nworkers:end] for i in 1:nworkers]
		for i in 2:nworkers
			MPI.send(part[i], comm; dest = (i - 1))
		end
		todo = part[1]
		starttime()
		#println("rank $rank | has indices $(todo)")
	else
		todo = MPI.recv(comm) # recieve indices
		#println("rank $rank | has indices $(todo)")
	end
	MPI.Barrier(comm)
	########################

	if Threads.nthreads() > length(todo)
		println("rank $rank | Warning: more threads than work")
	end
	if Threads.nthreads() <= length(todo)
		println("rank $rank | running with $(Threads.nthreads()) threads | $(length(todo)) trajectories to sample...")
	end

	### do work ###
	Threads.@threads for this_work in todo
		trajectory_id, circuit_id = this_work
		circuit = circuits[circuit_id]
		filename = "cluster/$(outputname)/idx$(circuit_id)_$(trajectory_id).jld2"
		sample_full(circuit, thermalization_time, n_samples, sample_distance; filename = filename)
		println("rank $rank | trajectory $(trajectory_id) of circuit $(circuit_id) -- done")
	end
	###############

	MPI.Barrier(comm)

	if rank == root
		folder = "cluster/$(outputname)"
		files = readdir(folder)
		to_read = filter(contains("idx"), files)
		idx_array = sort(unique(parse.(Int64, [split(split(string, "idx")[2], "_")[1] for string in to_read])))

		Threads.@threads for idx in idx_array
			this_files = filter(contains("idx$(idx)_"), to_read)
			I3s = zeros(length(this_files))
			EEs = zeros(length(this_files), length(subsystem_labels(circuits[1])))
			for i in eachindex(this_files)
				file = this_files[i]
				I3s[i] = jldopen("$(folder)/$(file)")["I3"]
				EEs[i, :] = jldopen("$(folder)/$(file)")["EE"]
				rm("$(folder)/$(file)")
			end
			writedlm("$(folder)/idx$(idx)_I3.txt", I3s)
			writedlm("$(folder)/idx$(idx)_EE.txt", EEs)
		end
		finaltime()
		MPI.Finalize()
	end
end

function mpi_sample_purification(
	circuits::Vector{T} where T <: AbstractCircuit,
	n_trajectories::Int64,
	timesteps::Int64,
	outputname::String,
)

	### MPI startup ###
	MPI.Init()
	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	nworkers = MPI.Comm_size(comm)
	root = 0
	MPI.Barrier(comm)
	###################

	### distribute work ###
	if rank == root
		meta = Dict(
			"circuits" => circuits,
			"n_trajectories" => n_trajectories,
			"timesteps" => timesteps,
			"outputname" => outputname,
		)
		meta_init(meta, outputname)
		work = distribute_work(circuits, n_trajectories)
		part = [work[i:nworkers:end] for i in 1:nworkers]
		for i in 2:nworkers
			MPI.send(part[i], comm; dest = (i - 1))
		end
		todo = part[1]
		starttime()
		#println("rank $rank | has indices $(todo)")
	else
		todo = MPI.recv(comm) # recieve indices
		#println("rank $rank | has indices $(todo)")
	end
	MPI.Barrier(comm)
	########################

	if Threads.nthreads() > length(todo)
		println("rank $rank | Warning: more threads than work")
	end
	if Threads.nthreads() <= length(todo)
		println("rank $rank | running with $(Threads.nthreads()) threads | $(length(todo)) trajectories to sample...")
	end

	### do work ###
	Threads.@threads for this_work in todo
		trajectory_id, circuit_id = this_work
		circuit = circuits[circuit_id]
		filename = "cluster/$(outputname)/idx$(circuit_id)_$(trajectory_id).jld2"
		sample_purification(circuit, timesteps; filename = filename)
		println("rank $rank | trajectory $(trajectory_id) of circuit $(circuit_id) -- done")
	end
	###############

	MPI.Barrier(comm)

	if rank == root
		folder = "cluster/$(outputname)"
		sample_purification_cleanup(folder, timesteps)
		finaltime()
		MPI.Finalize()
	end
end
