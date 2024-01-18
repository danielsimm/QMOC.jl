using QMOC
using QuantumClifford
using BenchmarkTools

L = 24



stabilizer = Stabilizer(inits)
destabilizer = Destabilizer(stabilizer)
mixeddestabilizer = MixedDestabilizer(stabilizer)

operators = QMOC.get_operators(QMOC.YaoKivelsonOrientableTest(L))
traj = QMOC.YaoKivelsonOrientableTest(L)
@benchmark QMOC.circuit!(stabilizer, $traj, $operators)
@benchmark QMOC.circuit!(destabilizer, $traj, $operators)
@benchmark QMOC.circuit!(mixeddestabilizer, $traj, $operators)

stab_tmi_time = @btime QMOC.tmi(stabilizer, $traj)
@benchmark QMOC.tmi(destabilizer, $traj)
@benchmark QMOC.tmi(mixeddestabilizer, $traj)

