using QMOC

L = 24
Kc = 0.654
circuits = [QMOC.YaoKivelsonNonorientableCircuit(L, [0.2, 0.8]), QMOC.YaoKivelsonNonorientableCircuit(L, [1-Kc, Kc]), QMOC.YaoKivelsonNonorientableCircuit(L, [0.5, 0.5])]

n_trajectories = 3840
timesteps = 3000


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"YKFullPurification24"
) 