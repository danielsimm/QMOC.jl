using QMOC
parameter_set = [[1/(1+K), K/(1+K)] for K in 0:0.5:10]
sim = simulation(:YaoKivelsonNonorientable, 12, "YaoKivelsonNonorientable_12", 200, parameter_set; checkpoints=false, verbosity=:info)
simulate(sim)