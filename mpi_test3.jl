using QMOC

L = 24

circuits = [QMOC.KitaevCircuit(L, [0.95, 0.025, 0.025])]

n_trajectories = 100
timesteps = L^2


QMOC.mpi_sample_purification(
	circuits,
	n_trajectories,
	timesteps,
	"purification_test2"
) 
