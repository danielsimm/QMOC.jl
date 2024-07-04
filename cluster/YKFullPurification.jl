using QMOC

L = 24
Kc = 0.6542
circuits = [QMOC.YaoKivelsonNonorientableCircuit(L, [0.2, 0.8]), QMOC.YaoKivelsonNonorientableCircuit(L, [1-Kc, Kc]), QMOC.YaoKivelsonNonorientableCircuit(L, [0.8, 0.2])]

n_trajectories = 3840
timesteps = 4000


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"YKFullPurification24"
) 