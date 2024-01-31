using QMOC

L = 12
function K_range(L, nu, Kc, ΔK)
    K_min = max(0, Kc - ΔK/(L^(1/nu)))
    K_max = Kc + ΔK/(L^(1/nu))
    return LinRange(K_min, K_max, 50)
end
Kc = 0.631
nu = 1.0
num_trajectories = 100
number_of_measurements = 100
Ks = K_range(L, nu, Kc, 3)
parameter_set = [[1-K, K] for K in Ks]
sim = simulation(:YaoKivelsonOrientable, L, "test", num_trajectories, parameter_set;  thermalization_steps=2*L, checkpoints=false, verbosity=:none, number_of_measurements=number_of_measurements)
simulate(sim)