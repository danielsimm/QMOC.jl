function _DHC_largeloop_sites(L)
    loops = []
    for i in 1:6*L^2
        if mod1(i, 6) == 3
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            site4 = _DHC_yneighbour(site3, L)
            site5 = _DHC_zneighbour(site4, L)
            site6 = _DHC_xneighbour(site5, L)
            site7 = _DHC_yneighbour(site6, L)
            site8 = _DHC_zneighbour(site7, L)
            site9 = _DHC_xneighbour(site8, L)
            site10 = _DHC_yneighbour(site9, L)
            site11 = _DHC_zneighbour(site10, L)
            site12 = _DHC_xneighbour(site11, L)
            push!(loops, [site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12])
        end
    end
    return loops
end

function _DHC_smallloop_sites(L)
    loops = []
    for i in 1:6*L^2
        if mod1(i, 6) == 2
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            push!(loops, [site1, site2, site3])
        elseif mod1(i, 6) == 4
            site1 = i
            site2 = _DHC_zneighbour(site1, L)
            site3 = _DHC_xneighbour(site2, L)
            push!(loops, [site1, site2, site3])
        end
    end
    return loops
end

function _DHC_wilsonline_sites(L)
    loop1 = zeros(Int, 6*L)
    loop1[1] = 1
    for i in 2:6*L
        if iseven(i)
            loop1[i] = _DHC_xneighbour(loop1[i-1], L)
        else
            loop1[i] = _DHC_zneighbour(loop1[i-1], L)
        end
    end

    loop2 = zeros(Int, 6*L)
    loop2[1] = 1
    for i in 2:6*L
        if iseven(i)
            loop2[i] = _DHC_xneighbour(loop2[i-1], L)
        else
            loop2[i] = _DHC_yneighbour(loop2[i-1], L)
        end
    end

    return loop1, loop2
end

function _DHC_largeloop_operators(L) :: Vector{PauliOperator}
    loops = _DHC_largeloop_sites(L)
    N = 6*L^2
    operators = []
    for loop in loops
        Xarr = falses(N)
        Zarr = falses(N)

        Xarr[loop[1]] = true

        Xarr[loop[2]] = true
        Zarr[loop[2]] = true

        Zarr[loop[3]] = true

        Xarr[loop[4]] = true

        Xarr[loop[5]] = true
        Zarr[loop[5]] = true

        Zarr[loop[6]] = true

        Xarr[loop[7]] = true

        Xarr[loop[8]] = true
        Zarr[loop[8]] = true

        Zarr[loop[9]] = true

        Xarr[loop[10]] = true

        Xarr[loop[11]] = true
        Zarr[loop[11]] = true

        Zarr[loop[12]] = true
       push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_smallloop_operators(L) :: Vector{PauliOperator}
    loops = _DHC_smallloop_sites(L)
    N = 6*L^2
    operators = []
    for loop in loops
        Xarr = falses(N)
        Zarr = falses(N)

        Xarr[loop[1]] = true

        Xarr[loop[2]] = true
        Zarr[loop[2]] = true

        Zarr[loop[3]] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_wilsonline_operators(L) :: Vector{PauliOperator}
    loop1, loop2 = _DHC_wilsonline_sites(L)
    N = 6*L^2
    operators = Vector{PauliOperator}(undef, 2)
    Xarr = falses(N)
    Zarr = falses(N)
    for i in eachindex(loop1)
        site = loop1[i]
        Xarr[site] = true
        Zarr[site] = true
    end
    operators[1] = PauliOperator(0x00, Xarr, Zarr)
    Xarr = falses(N)
    Zarr = falses(N)
    for i in eachindex(loop2)
        site = loop2[i]
        Zarr[site] = true
    end
    operators[2] = PauliOperator(0x00, Xarr, Zarr)
    return operators
end

function _DHC_ZZ_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:2:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Zarr[i] = true
        Zarr[_DHC_zneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return operators
end

function _DHC_XX_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Xarr[i] = true
        Xarr[_DHC_xneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return unique(operators)
end

function _DHC_YY_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    operators = []
    for i in 1:N # 3:6:N
        Xarr = falses(N)
        Zarr = falses(N)
        Xarr[i] = true
        Xarr[_DHC_yneighbour(i,L)] = true
        Zarr[i] = true
        Zarr[_DHC_yneighbour(i,L)] = true
        push!(operators, PauliOperator(0x00, Xarr, Zarr))
    end
    return unique(operators)
end

function _DHC_J_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    loops = _DHC_largeloop_sites(L)
    operators = []
    for loop in loops
        for direction in [:X, :Y, :Z]
            Xarr = falses(N)
            Zarr = falses(N)
            if direction == :X
                Xarr[loop[2]] = true
                Xarr[loop[3]] = true
            elseif direction == :Y
                Xarr[loop[1]] = true
                Xarr[loop[3]] = true
                Zarr[loop[1]] = true
                Zarr[loop[3]] = true
            elseif direction == :Z
                Zarr[loop[1]] = true
                Zarr[loop[2]] = true
            end
            push!(operators, PauliOperator(0x00, Xarr, Zarr))
        end
    end
    return unique(operators)
end

function _DHC_Jprime_operators(L) :: Vector{PauliOperator}
    N = 6*L^2
    loops = _DHC_largeloop_sites(L)
    operators = []
    for loop in loops
        for direction in [:X1, :X2, :Y1, :Y2, :Z1, :Z2]
            Xarr = falses(N)
            Zarr = falses(N)
            if direction == :Z1
                Zarr[loop[1]] = true
                Zarr[loop[2]] = true
            elseif direction == :Y1
                Xarr[loop[3]] = true
                Zarr[loop[3]] = true
                Xarr[loop[4]] = true
                Zarr[loop[4]] = true
            elseif direction == :X1
                Xarr[loop[5]] = true
                Xarr[loop[6]] = true
            elseif direction == :Z2
                Zarr[loop[7]] = true
                Zarr[loop[8]] = true
            elseif direction == :Y2
                Xarr[loop[9]] = true
                Zarr[loop[9]] = true
                Xarr[loop[10]] = true
                Zarr[loop[10]] = true
            elseif direction == :X2
                Xarr[loop[11]] = true
                Xarr[loop[12]] = true
            end
            push!(operators, PauliOperator(0x00, Xarr, Zarr))
        end
    end
    return unique(operators)
end