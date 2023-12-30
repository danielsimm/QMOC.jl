using QMOC
using QuantumClifford

function entropy(state::QuantumClifford.MixedDestabilizer, trajectory::DecoratedHoneycombTrajectory)
	algo = Val(:rref)
	L = trajectory.size
	EE = zeros(L + 1)
	for i in 1:L
        subsystem_range = 1:(i*6*L)
		EE[i+1] = entanglement_entropy(state, QMOC.DHC_subsystem(), algo)
	end
	return EE
end

function entropy_clipped(state::QuantumClifford.MixedDestabilizer, trajectory::DecoratedHoneycombTrajectory)
	bg = bigram(state; clip = true)
	L = trajectory.size
	EE = zeros(L + 1)
    for i in 1:L
        subsystem_range = 1:(i*6*L)
        EE[i+1] = length(subsystem_range) - count(r->(r[1] in subsystem_range && r[2] in subsystem_range), eachrow(bg))
    end
    return EE
end

p = 0.5
trajectory = QMOC.trajectory(
    :YaoKivelsonNonorientable, 
    12, 
    QMOC._number_of_qubits(:YaoKivelsonNonorientable, 2),
    "YaoKivelsonNonorientableTest",
    [p,1-p],
    false,
    :debug,
    1,
    3 * 12,
    3,
    1
)

state, time, existing_measurements = QMOC.hot_start(trajectory)

operators = QMOC.get_operators(trajectory)

state, time = QMOC.thermalise(state, trajectory, time, operators)

using BenchmarkTools
@benchmark entropy(state, trajectory)
@benchmark entropy_clipped(state, trajectory)
@benchmark QMOC.entropy(state, trajectory)
@benchmark QMOC.tmi(state, trajectory)