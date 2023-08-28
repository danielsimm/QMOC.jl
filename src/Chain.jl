function chain_timestep!(stabilizer, px, py, pz, type) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    if type == :PP
       _chain_PP!(stabilizer, px, py, pz)
    elseif type == :PQ
        _chain_PQ!(stabilizer, px, py, pz)
    else
        error("type must be :PP or :PQ")
    end
end


function _chain_PP!(stabilizer, px, py, pz) # measures XX, YY, ZZ on neighbouring sites (first site chosen randomly) with probability px, py, pz respectively
    numberOfQubits = nqubits(stabilizer)

    for subtime in 1:numberOfQubits
        firstsite = rand(1:numberOfQubits)
        secondsite = mod1(firstsite+1, numberOfQubits)
        
        X_arr = falses(numberOfQubits)
        Z_arr = falses(numberOfQubits)

        probability = rand()

        if probability < px
            X_arr[firstsite] = true
            X_arr[secondsite] = true
        elseif probability < (px + py)
            X_arr[firstsite] = true
            X_arr[secondsite] = true
            Z_arr[firstsite] = true
            Z_arr[secondsite] = true
        else
            Z_arr[firstsite] = true
            Z_arr[secondsite] = true
        end

        project!(stabilizer, PauliOperator(0x00, X_arr,Z_arr), phases=false, keep_result=false)
    end
end

function _chain_PQ!(stabilizer, px, py, pz)
    numberOfQubits = nqubits(stabilizer)

    for subtime in 1:numberOfQubits
        site1 = rand(1:numberOfQubits)
        site2 = mod1(site1+1, numberOfQubits)
        X_arr = falses(numberOfQubits)
        Z_arr = falses(numberOfQubits)

        probability1 = rand()
        if probability1 < px
            X_arr[site1] = true
        elseif probability1 < (px + py)
            X_arr[site1] = true
            Z_arr[site1] = true
        else
            Z_arr[site1] = true
        end
        probability2 = rand()
        if probability2 < px
            X_arr[site2] = true
        elseif probability2 < (px + py)
            X_arr[site2] = true
            Z_arr[site2] = true
        else
            Z_arr[site2] = true
        end
        project!(stabilizer, PauliOperator(0x00, X_arr,Z_arr), phases=false, keep_result=false)
    end
end


### Observables

function Chain_EE(state; algo=Val(:rref))
    N = nqubits(state)
    EE = zeros(33)
    subsystems = Int.(collect(range(0, N, 33)))
    for i in 2:32
        EE[i] = entanglement_entropy(state, 1:subsystems[i], algo)
    end
    return EE
end

function Chain_TMI(state; algo=Val(:rref)) # no geometry
    N = nqubits(state)
    A = 1:Int(N/4)
    B = Int(N/4)+1:Int(N/2)
    C = Int(N/2)+1:Int(3N/4)
    
    if mod(N, 4) != 0
        @info "L must be a multiple of 4, but is $(N). Tripartite mutual information is ill-defined."
    end

    SA = entanglement_entropy(state, A, algo)
    SB = entanglement_entropy(state, B, algo)
    SC = entanglement_entropy(state, C, algo)
    SAB = entanglement_entropy(state, union(A,B), algo)
    SBC = entanglement_entropy(state, union(B,C), algo)
    SAC = entanglement_entropy(state, union(A,C), algo)
    SABC = entanglement_entropy(state, union(A, B, C), algo)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end