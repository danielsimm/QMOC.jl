using QMOC

L = 64
pc = 0.652
circuits = [QMOC.KitaevCircuit(L, [1/3, 1/3, 1/3]), QMOC.KitaevCircuit(L, [1/2, 1/2, 0.0]), QMOC.KitaevCircuit(L, [pc, (1-pc)/2, (1-pc)/2]), QMOC.KitaevCircuit(L, [0.7, 0.15, 0.15])]

n_trajectories = 3840
timesteps = 3*L


QMOC.mpi_sample_purification(
	circuits,
	n_trajectories,
	timesteps,
	"KitaevPurification64"
) 