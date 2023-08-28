function graphEE(graph, adjmat, subsystem)
    other_subsystem = filter(i->!(i in collect(subsystem)), 1:Graphs.nv(graph))
    subadjmat = Nemo.matrix(Nemo.residue_ring(Nemo.ZZ, 2), collect(adjmat[subsystem,other_subsystem]))
    return LinearAlgebra.rank(subadjmat)
end

# slow...
function observables(state::QuantumClifford.AbstractStabilizer, geometry::Symbol)
    N = nqubits(state)
    graph = Graphs.Graph(state)
    adjmat = Graphs.adjacency_matrix(graph)
    if geometry == :Chain
        ## EE sweep
        EE = zeros(33)
        subsystems = Int.(collect(range(0, N, 33)))
        for i in 2:32
            EE[i] = graphEE(graph, adjmat, 1:subsystems[i])
        end
        ## TMI
        A = 1:Int(N/4)
        B = Int(N/4)+1:Int(N/2)
        C = Int(N/2)+1:Int(3N/4)
        SA = graphEE(graph, adjmat, A)
        SB = graphEE(graph, adjmat, B)
        SC = graphEE(graph, adjmat, C)
        SAB = graphEE(graph, adjmat, union(A,B))
        SBC = graphEE(graph, adjmat, union(B,C))
        SAC = graphEE(graph, adjmat, union(A,C))
        SABC = graphEE(graph, adjmat, union(A, B, C))
        TMI = SA + SB + SC - SAB - SBC - SAC + SABC
    elseif geometry == :Honeycomb
        L = Int(sqrt(N/2))
        subsystem = HC_subsystem
        EE = zeros(L+1)
    elseif geometry == :DecoratedHoneycomb
        L = Int(sqrt(N/3))
        subsystem = DHC_subsystem
        EE = zeros(L+1)
    else
        error("Geometry not supported")
    end

    if geometry != :Chain
        for i in 1:L
            EE[i+1] = graphEE(graph, adjmat, subsystem(L, 1:i))
        end
        A = subsystem(L, 1:Int(L/4))
        B = subsystem(L, Int(L/4)+1:Int(L/2))
        C = subsystem(L, Int(L/2)+1:Int(3L/4))
        AB = union(A,B)
        BC = union(B,C)
        AC = union(A,C)
        ABC = union(A,B,C)
        SA = graphEE(graph, adjmat, A)
        SB = graphEE(graph, adjmat, B)
        SC = graphEE(graph, adjmat, C)
        SAB = graphEE(graph, adjmat, AB)
        SBC = graphEE(graph, adjmat, BC)
        SAC = graphEE(graph, adjmat, AC)
        SABC = graphEE(graph, adjmat, ABC)
        TMI = SA + SB + SC - SAB - SBC - SAC + SABC
    end
    return EE, TMI
end