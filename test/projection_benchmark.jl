using QMOC
using QuantumClifford
using BenchmarkTools

function _DHC_J_sites_nonorientable(L)
    length = 2*L^2
    loops = QMOC._DHC_smallloop_sites(L)
    X_sites = [zeros(Int, 2) for _ in 1:length]
    Y_sites = [zeros(Int, 2) for _ in 1:length]
    Z_sites = [zeros(Int, 2) for _ in 1:length]
    for (i, loop) in enumerate(loops)
        for direction in [:X, :Y, :Z]
            if direction == :X
                X_sites[i] = [loop[2], loop[3]]
            elseif direction == :Y
                Y_sites[i] = [loop[1], loop[3]]
            elseif direction == :Z
                Z_sites[i] = [loop[1], loop[2]]
            end
        end
    end
    return X_sites, Y_sites, Z_sites
end

function _DHC_J_sites_orientable(L)
    length = 2*L^2
    loops = QMOC._DHC_smallloop_sites(L)
    X_sites = [zeros(Int, 2) for _ in 1:length]
    Y_sites = [zeros(Int, 2) for _ in 1:length]
    for (i, loop) in enumerate(loops)
        for direction in [:X, :Y]
            if direction == :X
                X_sites[i] = [loop[2], loop[3]]
            elseif direction == :Y
                Y_sites[i] = [loop[1], loop[3]]
            end
        end
    end
    return X_sites, Y_sites
end

function _DHC_K_sites(L) # every bond exists twice, in both directions - irrelevant for sampling
    length = 2*L^2
    offset = L^2
    loops = QMOC._DHC_largeloop_sites(L)
    X_sites = [zeros(Int, 2) for _ in 1:length]
    Y_sites = [zeros(Int, 2) for _ in 1:length]
    Z_sites = [zeros(Int, 2) for _ in 1:length]
    for (i, loop) in enumerate(loops)
        for direction in [:X1, :X2, :Y1, :Y2, :Z1, :Z2]
            if direction == :Z1
                Z_sites[i] = [loop[1], loop[2]]
            elseif direction == :Y1
                Y_sites[i] = [loop[3], loop[4]]
                Y_sites[i] = [loop[3], loop[4]]
            elseif direction == :X1
                X_sites[i] = [loop[5], loop[6]]
            elseif direction == :Z2
                Z_sites[i+offset] = [loop[7], loop[8]]
            elseif direction == :Y2
                Y_sites[i+offset] = [loop[9], loop[10]]
                Y_sites[i+offset] = [loop[9], loop[10]]
            elseif direction == :X2
                X_sites[i+offset] = [loop[11], loop[12]]
            end
        end
    end
    return X_sites, Y_sites, Z_sites
end


L = 48
traj = QMOC.YaoKivelsonOrientableTest(L)
K_sites_X, K_sites_Y, K_sites_Z = _DHC_K_sites(L)
J_sites_X, J_sites_Y = _DHC_J_sites_orientable(L)
inits = QMOC.initialise(traj)

stabilizer = Stabilizer(inits)
destabilizer = Destabilizer(stabilizer)
mixeddestabilizer = MixedDestabilizer(stabilizer)

function alternative_circuit!(state::QuantumClifford.MixedDestabilizer, trajectory::YaoKivelsonOrientableTrajectory, J_sites_X, J_sites_Y, K_sites_X, K_sites_Y, K_sites_Z) ::Nothing
    J = trajectory.params[1] # probability of triangular measurement, orientable: silence triangular Z bond, non-orientable: all bonds
    # K = trajectory.params[2] # probability of hexagonal/Kitaev-style measurement
    L = trajectory.size
    J_length = 2*L^2
    K_length = 2*L^2
    for subtime in 1:trajectory.nqubits
        if rand() < J # pick random J site and direction
            dir = rand([:X, :Y])
            if dir == :X
                site = J_sites_X[rand(1:J_length)]
                apply!(state, sCNOT(site[1], site[2]))
                projectX!(state, site[1], keep_result=false, phases=false)
                apply!(state, sCNOT(site[1], site[2]))
            else
                site = J_sites_Y[rand(1:J_length)]
                apply!(state, sCNOT(site[1], site[2]))
                projectY!(state, site[1], keep_result=false, phases=false)
                apply!(state, sCNOT(site[1], site[2]))
            end
        else # pick random K site and direction
            dir = rand([:X, :Y, :Z])
            if dir == :X
                site = K_sites_X[rand(1:K_length)]
                apply!(state, sCNOT(site[1], site[2]))
                projectX!(state, site[1], keep_result=false, phases=false)
                apply!(state, sCNOT(site[1], site[2]))
            elseif dir == :Y
                site = K_sites_Y[rand(1:K_length)]
                apply!(state, sCNOT(site[1], site[2]))
                projectY!(state, site[1], keep_result=false, phases=false)
                apply!(state, sCNOT(site[1], site[2]))
            else
                site = K_sites_Z[rand(1:K_length)]
                apply!(state, sCNOT(site[1], site[2]))
                projectZ!(state, site[1], keep_result=false, phases=false)
                apply!(state, sCNOT(site[1], site[2]))
            end
        end
    end
    return nothing
end
operators = QMOC.get_operators(QMOC.YaoKivelsonOrientableTest(L))

@benchmark QMOC.circuit!(stabilizer, $traj, $operators)
@benchmark QMOC.circuit!(destabilizer, $traj, $operators)
@benchmark QMOC.circuit!(mixeddestabilizer, $traj, $operators)
@benchmark alternative_circuit!(mixeddestabilizer, $traj, $J_sites_X, $J_sites_Y, $K_sites_X, $K_sites_Y, $K_sites_Z)

stab_tmi_time = @btime QMOC.tmi(stabilizer, $traj)
@benchmark QMOC.tmi(destabilizer, $traj)
@benchmark QMOC.tmi(mixeddestabilizer, $traj)

