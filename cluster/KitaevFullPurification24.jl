using QMOC

L = 24
pc = 0.6524
circuits = [QMOC.KitaevCircuit(L, [1/3, 1/3, 1/3])], QMOC.KitaevCircuit(L, [1/2, 1/2, 0.0]), QMOC.KitaevCircuit(L, [pc, (1-pc)/2, (1-pc)/2]), QMOC.KitaevCircuit(L, [0.7, 0.15, 0.15])]

n_trajectories = 3840
timesteps = 2*L^2*1000


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"KitaevFullPurificationDetail24"
) 