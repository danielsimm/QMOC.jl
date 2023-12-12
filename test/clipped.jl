using QMOC
traj = QMOC.trajectory(:YaoKivelsonOrientable, 24, QMOC._number_of_qubits(:YaoKivelsonOrientable, 24), "test", [1//3, 1//3, 1//3], true, :debug, 1, 3*24, 3, 100)
state = QMOC.initialise(traj)
operators = QMOC.get_operators(traj)
t = 0
state, t = QMOC.thermalise(state, traj, t, operators)

QMOC.DHC_subsystem(24, [1])

