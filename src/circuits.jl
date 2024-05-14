import Base: hash
import Graphs
import SimpleWeightedGraphs as swg

abstract type AbstractCircuit end
abstract type DecoratedHoneycombCircuit <: AbstractCircuit end
abstract type HoneycombCircuit <: AbstractCircuit end
abstract type ChainCircuit <: AbstractCircuit end

"""
    hash(c::AbstractCircuit) -> UInt

    Returns a hash of the trajectory properties.
"""
function hash(c::AbstractCircuit)
    return hash("$(typeof(c))_$(c.size)_$(c.params)")
end
struct GenericCircuit <: AbstractCircuit
    size::Int
    dims::Int
    nqubits::Int
    graph::swg.SimpleWeightedGraph
    name::String
    function GenericCircuit(size, dims, nqubits, graph, name)
        return new(size, dims, nqubits, graph, name)
    end
end

struct YaoKivelsonXYZCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonXYZCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct YaoKivelsonOrientableCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonOrientableCircuit(size, J::Real)
        return new(size, 6*size^2, [J, 1-J])
    end
    function YaoKivelsonOrientableCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct YaoKivelsonNonorientableCircuit <: DecoratedHoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function YaoKivelsonNonorientableCircuit(size, J::Real)
        return new(size, 6*size^2, [J, 1-J])
    end
    function YaoKivelsonNonorientableCircuit(size, params)
        return new(size, 6*size^2, params)
    end
end

struct KitaevCircuit <: HoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function KitaevCircuit(size, params)
        return new(size, 2*size^2, params)
    end
end

struct KekuleCircuit <: HoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function KekuleCircuit(size, params)
        if length(params) == 3
            return new(size, 2*size^2, params)
        else
            error("KekuleCircuit requires 3 parameters Jˣ, Jʸ, Jᶻ")
        end
    end
end

struct KitaevNNNCircuit <: HoneycombCircuit
    size::Int
    nqubits::Int
    params::Vector{Real}
    function KitaevNNNCircuit(size, params)
        if length(params) == 6
            if sum(params) ≈ 1
                return new(size, 2*size^2, params)
            else
                error("KitaevNNNCircuit requires parameters to sum to 1")
            end
        else
            error("KitaevNNNCircuit requires 6 parameters Kˣ, Kʸ, Kᶻ, Jˣ, Jʸ, Jᶻ")
        end
    end
end

struct AncillaCircuit <: AbstractCircuit
    circuit::AbstractCircuit
    ancillas::Vector{Int}
end

function _mutual_information(stab, ind1, ind2)
    SA = entanglement_entropy(stab, [ind1], Val(:rref))
    SB = entanglement_entropy(stab, [ind2], Val(:rref))
    SAB = entanglement_entropy(stab, [ind1, ind2], Val(:rref))
    return SA + SB - SAB
end

function _mutual_information(stab, A::Vector{Int}, B::Vector{Int})
    SA = entanglement_entropy(stab, A, Val(:rref))
    SB = entanglement_entropy(stab, B, Val(:rref))
    SAB = entanglement_entropy(stab, union(A, B), Val(:rref))
    return SA + SB - SAB
end


function DHC_rows(L)
    rows = [zeros(Int, 6*L) for _ in 1:L]
    for i in eachindex(rows)
       rows[i][1] = (i-1)*6*L+1
       for j in 2:6*L
        if iseven(j)
            rows[i][j] = _DHC_xneighbour(rows[i][j-1], L)
        else
            rows[i][j] = _DHC_yneighbour(rows[i][j-1], L)
        end
       end
    end
    return rows
end

function DHC_ancilla_sites(L)
    rows = DHC_rows(L)
    ancilla1 = 1
    # row 1+L/2 col 1+L/2
    shift = L-div(L,2)
    ancilla2 = rows[1+div(L,2)][3*L+1]
    return [ancilla1, ancilla2]
end

DHC_snake(L) = reduce(vcat, DHC_rows(L))


function DHC_ancilla_strings(L, systemqubits)
    sites = DHC_snake(L)
    N = 6*L^2
    stringqubits1 = sites[findfirst(isequal(systemqubits[1]), sites)+1:end]
    stringqubits2 = sites[findfirst(isequal(systemqubits[2]), sites)+1:end]
    Zarr1 = zeros(Bool, N + 2)
    Zarr2 = zeros(Bool, N + 2)
    Xarr1 = zeros(Bool, N + 2)
    Xarr2 = zeros(Bool, N + 2)
    for ind in stringqubits1
        Zarr1[ind] = true
    end
    for ind in stringqubits2
        Zarr2[ind] = true
    end
    Xarr1[systemqubits[1]] = true
    Xarr1[N+1] = true
    Xarr2[systemqubits[2]] = true
    Xarr2[N+2] = true
    Zarr1b = deepcopy(Zarr1)
    Zarr2b = deepcopy(Zarr2)
    Zarr1b[systemqubits[1]] = true
    Zarr1b[N+1] = true
    Zarr2b[systemqubits[2]] = true
    Zarr2b[N+2] = true
    # PauliOperator(zeros(Bool, N+2), Zarr1b), PauliOperator(zeros(Bool, N+2), Zarr2b),
    return [PauliOperator(Xarr1, Zarr1), PauliOperator(Xarr2, Zarr2)]
end

function add_ancillas(stabilizer, circuit)
    L = circuit.size
    N = circuit.nqubits
    new_stabilizer_paulis = [PauliOperator(zeros(Bool, N+2), zeros(Bool, N+2)) for _ in 1:N+2]
    for i in 1:N
        new_stabilizer_paulis[i] = embed(N+2, 1:N, stabilizerview(stabilizer)[i])
    end
    ancilla1Zarr = zeros(Bool, N+2)
    ancilla1Xarr = zeros(Bool, N+2)
    ancilla2Zarr = zeros(Bool, N+2)
    ancilla2Xarr = zeros(Bool, N+2)
    ancilla1Zarr[N+1] = true
    ancilla2Zarr[N+2] = true
    new_stabilizer_paulis[N+1] = PauliOperator(ancilla1Xarr, ancilla1Zarr)
    new_stabilizer_paulis[N+2] = PauliOperator(ancilla2Xarr, ancilla2Zarr)

    return Stabilizer(new_stabilizer_paulis)
end

function ancilla_scheme(circuit, n_subsamples=100)
    L = circuit.size
    N = circuit.nqubits
    stab = thermal_state(circuit)

    # add ancillas
    stab = add_ancillas(stab, circuit)
    coupling_ops = DHC_ancilla_strings(L, DHC_ancilla_sites(L))
    for op in coupling_ops
        project!(stab, op, keep_result=false, phases=false)
    end


    evolution_ops = get_operators(circuit)
    for i in eachindex(evolution_ops)
        evolution_ops[i] = QuantumClifford.embed(N+2, 1:N, evolution_ops[i])
    end

    # evolve system
    for t in 1:2*L
        apply!(stab, circuit, evolution_ops)
    end

    subsamples = zeros(Int, n_subsamples)
    for subsample in 1:n_subsamples
        for _ in 1:3
            apply!(stab, circuit, evolution_ops)
        end
        eval_stab = deepcopy(stab)
        for op in _DHC_K_operators(L)
            this_op = QuantumClifford.embed(N+2, 1:N, op)
            project!(eval_stab, this_op, keep_result=false, phases=false)
        end
        subsamples[subsample] = _mutual_information(eval_stab, N+1, N+2)
    end

    return sum(subsamples)/n_subsamples

end

