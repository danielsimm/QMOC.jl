using QMOC

L = 24
Kc = 0.654
circuits = [QMOC.YaoKivelsonNonorientableCircuit(L, Kc-0.25), QMOC.YaoKivelsonNonorientableCircuit(L, Kc), QMOC.YaoKivelsonNonorientableCircuit(L, Kc+0.25)]

n_trajectories = 3840
timesteps = L^3


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"YKFullPurification24"
) 