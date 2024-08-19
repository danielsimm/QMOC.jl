using QMOC

L = 24
Kc = 0.6542
circuits = [QMOC.YaoKivelsonNonorientableCircuit(L, [0.25, 0.75]), QMOC.YaoKivelsonNonorientableCircuit(L, [1-Kc, Kc]), QMOC.YaoKivelsonNonorientableCircuit(L, [0.45, 0.55])]

n_trajectories = 3840
timesteps = 6*L^2*150


QMOC.mpi_sample_full_purification(
	circuits,
	n_trajectories,
	timesteps,
	"YKFullPurificationDetail24"
) 