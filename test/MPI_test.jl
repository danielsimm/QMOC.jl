using QMOC
parameter_set = [[1/(1+K), K/(1+K)] for K in 1.4:0.05:2.4]
parameter_set2 = [[1/(1+K), K/(1+K)] for K in 0.0:0.1:5.0]
parameter_set = union(parameter_set, parameter_set2)
sim = simulation(:YaoKivelsonNonorientable, 12, "test", 10, parameter_set2; checkpoints=false, verbosity=:none)
simulate(sim)