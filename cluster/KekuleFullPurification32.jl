using QMOC

L = 32
pc = 0.623
circuits = [QMOC.KekuleCircuit(L, [1/3, 1/3, 1/3]), QMOC.KekuleCircuit(L, [pc, (1-pc)/2, (1-pc)/2]), QMOC.KekuleCircuit(L, [0.7, 0.15, 0.15])]

n_trajectories = 3840
timesteps = 3000


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"KekuleFullPurification32"
) 