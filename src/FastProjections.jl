function projectXX!(state::QuantumClifford.AbstractStabilizer, site1::Int, site2::Int)
    apply!(state, sCNOT(site1, site2))
    projectX!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end

function projectYY!(state::QuantumClifford.AbstractStabilizer, site1::Int, site2::Int)
    apply!(state, sCNOT(site1, site2))
    projectY!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end

function projectZZ!(state::QuantumClifford.AbstractStabilizer, site1::Int, site2::Int)
    apply!(state, sCNOT(site1, site2))
    projectZ!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end

function randomXX!(state::MixedDestabilizer, trajectory::PPChainTrajectoryFast)
    site1 = rand(1:trajectory.nqubits)
    site2 = mod1(site1+1, trajectory.nqubits)
    apply!(state, sCNOT(site1, site2))
    projectX!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end

function randomYY!(state::MixedDestabilizer, trajectory::PPChainTrajectoryFast)
    site1 = rand(1:trajectory.nqubits)
    site2 = mod1(site1+1, trajectory.nqubits)
    apply!(state, sCNOT(site1, site2))
    projectY!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end

function randomZZ!(state::MixedDestabilizer, trajectory::PPChainTrajectoryFast)
    site1 = rand(1:trajectory.nqubits)
    site2 = mod1(site1+1, trajectory.nqubits)
    apply!(state, sCNOT(site1, site2))
    projectZ!(state, site2, keep_result=false, phases=false)
    apply!(state, sCNOT(site1, site2))
end