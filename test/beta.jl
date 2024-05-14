using QMOC
import QuantumClifford as qc
L = 16
N = 6 * L^2
function K_range(L, nu, Kc, ΔK)
	K_min = max(0, Kc - ΔK / (L^(1 / nu)))
	K_max = Kc + ΔK / (L^(1 / nu))
	return LinRange(K_min, K_max, 16)
end
Kc_orientable = 0.631
nu_orientable = 1.0
Kc_nonorientable = 0.654
nu_nonorientable = 0.94
ΔK = 3
Ks_orientable = K_range(L, nu_orientable, Kc_orientable, ΔK)
Ks_nonorientable = K_range(L, nu_nonorientable, Kc_nonorientable, ΔK)
parameter_set_orientable = [[1 - K, K] for K in Ks_orientable]
parameter_set_nonorientable = [[1 - K, K] for K in Ks_nonorientable]

circuits = [QMOC.YaoKivelsonOrientableCircuit(L, p) for p in parameter_set_orientable]

mutual_information = zeros(length(circuits))

Threads.@threads for i in eachindex(circuits)
	c = circuits[i]
	this_mi = 0
	for _ in 1:100
		this_mi += QMOC.ancilla_scheme(c)
	end
	mutual_information[i] = this_mi / 30
end

for i in eachindex(circuits)
	println("$(parameter_set_orientable[i]) -> $(mutual_information[i])")
end

QMOC.DHC_ancilla_strings(4, QMOC.DHC_ancilla_sites(4))

test_state = QMOC.thermal_state(circuits[1])
test_state = QMOC.add_ancillas(test_state, circuits[1])
qc.stabilizerview(test_state)
qc.entanglement_entropy(test_state, [N + 1], Val(:rref))
N + 1


new_stabilizer_paulis = [qc.PauliOperator(zeros(Bool, N + 2), zeros(Bool, N + 2)) for _ in 1:N+2]
for i in 1:N
	new_stabilizer_paulis[i] = qc.embed(N + 2, 1:N, qc.stabilizerview(test_state)[i])
end
for i in N+1:N+2
	new_stabilizer_paulis[i] = qc.PauliOperator(zeros(Bool, N + 2), zeros(Bool, N + 2))
end
new_stabilizer_paulis

qc.Stabilizer(new_stabilizer_paulis) == qc.stabilizerview(qc.MixedDestabilizer(qc.Stabilizer(new_stabilizer_paulis)))
