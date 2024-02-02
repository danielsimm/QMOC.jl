using QMOC
using CairoMakie
### The initial state should fix all conserved qunatities in the +1 eigenstate, since they will inevitably be measured after infinitely long times
### The set of conserved quantities/stabilisers consists of the products of the operators on all dodecagons and all triangles
### Since we're generally working with periodic boundary conditions, we also need to include the products of the operators along two non-contractible loops
L = 8
largeloops = QMOC._DHC_largeloop_operators(L)
smallloops = QMOC._DHC_smallloop_operators(L)
wilsonlines = QMOC._DHC_wilsonline_operators(L)
bilinears = QMOC._DHC_ZZ_operators(L)
evolution_operators = [QMOC._DHC_J_operators_nonorientable(L)..., QMOC._DHC_K_operators(L)...]

### First, we confirm that the evolution operators commute with all conserved quantities and with each other
function commute_all(operatorsA, operatorsB)
    for operatorA in operatorsA
        for operatorB in operatorsB
            if QuantumClifford.comm(operatorA, operatorB) != 0
                return false
            end
        end
    end
    return true
end
commute_all(evolution_operators, largeloops)
commute_all(evolution_operators, smallloops)
commute_all(evolution_operators, wilsonlines)
commute_all(largeloops, smallloops)
commute_all(largeloops, wilsonlines)
commute_all(smallloops, wilsonlines)
commute_all(smallloops, smallloops)
commute_all(largeloops, largeloops)
commute_all(wilsonlines, wilsonlines)


function diagonal_init(L)
    return one(Stabilizer, 6*L^2)
end

function trivial_init(L)
    stab = one(Stabilizer, 6*L^2)
    largeloops = QMOC._DHC_largeloop_operators(L)
    smallloops = QMOC._DHC_smallloop_operators(L)
    wilsonlines = QMOC._DHC_wilsonline_operators(L)
    for op in [largeloops..., smallloops..., wilsonlines...]
        QuantumClifford.project!(stab, op, keep_result=false, phases=false)
    end
    return stab
end

function ZZ_init(L)
    stab = one(Stabilizer, 6*L^2)
    largeloops = QMOC._DHC_largeloop_operators(L)
    smallloops = QMOC._DHC_smallloop_operators(L)
    wilsonlines = QMOC._DHC_wilsonline_operators(L)
    bilinears = QMOC._DHC_ZZ_operators(L)
    for op in [largeloops..., smallloops..., wilsonlines..., bilinears...]
        QuantumClifford.project!(stab, op, keep_result=false, phases=false)
    end
    return stab
end

### The easiest way to construct the initial state is to start with a trivial stabilizer and then project all operators onto it
stab = trivial_init(L)
commute_all(wilsonlines, stab)
commute_all(largeloops, stab)
commute_all(smallloops, stab)

### We can also fix a configuration of the bilinears
stab = ZZ_init(L)
commute_all(wilsonlines, stab)
commute_all(largeloops, stab)
commute_all(smallloops, stab)

### We can test the time to reach the steady state by evolving the initial state with the evolution operators

function subsystem(rows, L)
    start = minimum(rows)
    stop = maximum(rows)
    return ((start-1)*6L)+1:stop*6L
end
function YK_tmi(state, L)
    algo=Val(:rref)
    if mod(L, 4) != 0
        @info "L must be a multiple of 4, but is $(L). Tripartite mutual information is ill-defined."
    end
    A = subsystem(1:Int(L/4),  L)
    B = subsystem(Int(L/4)+1:Int(L/2), L)
    C = subsystem(Int(L/2)+1:Int(3L/4), L)
    SA = entanglement_entropy(state, A, algo)
    SB = entanglement_entropy(state, B, algo)
    SC = entanglement_entropy(state, C, algo)
    SAB = entanglement_entropy(state, union(A,B), algo)
    SBC = entanglement_entropy(state, union(B,C), algo)
    SAC = entanglement_entropy(state, union(A,C), algo)
    SABC = entanglement_entropy(state, union(A, B, C), algo)
    return SA + SB + SC - SAB - SBC - SAC + SABC
end
function layer!(state, L, ops)
    for _ in 1:6*L^2
        project!(state, rand(ops), keep_result=false, phases=false)
    end
    return nothing
end
        
# diagonal_state
begin
    TMIs_diag = zeros(5*L+1)
    init = diagonal_init(L)
    for samples in 1:1000
        stab = copy(init)
        TMIs_diag[1] += YK_tmi(stab, L)
        for i in 1:5*L
            layer!(stab, L, evolution_operators)
            TMIs_diag[i+1] += YK_tmi(stab, L)
        end
    end
    TMIs_diag ./= 1000
end
begin
    TMIs_triv = zeros(5*L+1)
    init = trivial_init(L)
    for samples in 1:1000
        stab = copy(init)
        TMIs_triv[1] += YK_tmi(stab, L)
        for i in 1:5*L
            layer!(stab, L, evolution_operators)
            TMIs_triv[i+1] += YK_tmi(stab, L)
        end
    end
    TMIs_triv ./= 1000
end
with_theme(theme_latexfonts()) do 
    fig = Figure(resolution = (800, 600), font =:sans)
    ax = Axis(fig[1, 1], xlabel = "Layer", ylabel = "Tripartite mutual information", xticks = 0:5*L)
    lines!(ax, (0:5*L)./L, TMIs_triv, color = :black, linewidth = 2)
    fig
end

# ZZ_state
begin
    TMIs_ZZ = zeros(5*L+1)
    init = ZZ_init(L)
    for samples in 1:1000
        stab = copy(init)
        TMIs_ZZ[1] += YK_tmi(stab, L)
        for i in 1:5*L
            layer!(stab, L, evolution_operators)
            TMIs_ZZ[i+1] += YK_tmi(stab, L)
        end
    end
    TMIs_ZZ ./= 1000
end
with_theme(theme_latexfonts()) do 
    fig = Figure(resolution = (800, 600), font =:sans)
    ax = Axis(fig[1, 1], xlabel = "Layer/L", ylabel = "Tripartite mutual information", xticks = 0:5*L)
    lines!(ax, (0:5*L)./L, TMIs_ZZ, color = :blue, linewidth = 2, label = "ZZ")
    lines!(ax, (0:5*L)./L, TMIs_triv, color = :red, linewidth = 2, label = "Trivial")
    axislegend(ax, position = :rt)
    fig
end

init = trivial_init(16)
ops= [QMOC._DHC_J_operators_nonorientable(16)..., QMOC._DHC_K_operators(16)...]
entropy = zeros(17)
for sample in 1:100
    state = copy(init)
    for _ in 1:3*16
        layer!(state, 16, ops)
    end
    entropy[2:end] .+= [QuantumClifford.entanglement_entropy(state, subsystem(1:i, 16), Val(:rref)) for i in 1:16]
end

with_theme(theme_latexfonts()) do 
    fig = Figure(resolution = (800, 600), font =:sans)
    ax = Axis(fig[1,1])
    lines!(ax, 1:17, entropy./100, color = :black, linewidth = 2)
    fig
end

function DHC_snake(L)
    sites = zeros(Int, 6*L^2)
    rows = [((i-1)*6*L)+1:i*6*L for i in 1:L]
    for r in eachindex(rows)
        if iseven(r)
            sites[rows[r]] = rows[r][end:-1:1]
        else
            sites[rows[r]] = rows[r]
        end
    end
    return sites
end

DHC_snake(4)

function entropy_sweep(state, L)
    ranges = [subsystem(1:i, L) for i in 1:L]
    stab = copy(stabilizerview(state))
    QuantumClifford.permute!(stab, DHC_snake(L))
    bg = QuantumClifford.bigram(stab)
    [length(subsystem_range) - count(r->(r[1] in subsystem_range && r[2] in subsystem_range), eachrow(bg)) for subsystem_range in ranges]
end
using BenchmarkTools
L = 32
ops = [QMOC._DHC_J_operators_nonorientable(L)..., QMOC._DHC_K_operators(L)...]
init = trivial_init(L)
for _ in 1:3*L
    layer!(init, L, ops)
end
@benchmark entropy_sweep($init, 32)
@benchmark [QuantumClifford.entanglement_entropy($init, subsystem(1:i, 32), Val(:rref)) for i in 1:16]

