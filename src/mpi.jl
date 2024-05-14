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
		println("rank $rank | has indices $(todo)")
	else
		todo = MPI.recv(comm) # recieve indices
		println("rank $rank | has indices $(todo)")
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
	MPI.Finalize()
end