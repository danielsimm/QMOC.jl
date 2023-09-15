function _HC_ZZ_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j = _HC_zneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        Zarr[i] = true
        Zarr[j] = true
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
    return operators
end

function _HC_XX_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j = _HC_xneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        Xarr[i] = true
        Xarr[j] = true
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
    return operators
end

function _HC_YY_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j = _HC_yneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        Xarr[i] = true
        Xarr[j] = true
        Zarr[i] = true
        Zarr[j] = true
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
    return operators
end

function _HC_red_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j, direction = _HC_redneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        if direction == :X
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Y
            Zarr[i] = true
            Zarr[j] = true
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Z
            Zarr[i] = true
            Zarr[j] = true
        end
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
end

function _HC_green_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j, direction = _HC_greenneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        if direction == :X
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Y
            Zarr[i] = true
            Zarr[j] = true
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Z
            Zarr[i] = true
            Zarr[j] = true
        end
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
end

function _HC_blue_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2)
    for k in eachindex(operators)
        i = 2*k
        j, direction = _HC_blueneighbour(i, L)
        Xarr = falses(2*L^2)
        Zarr = falses(2*L^2)
        if direction == :X
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Y
            Zarr[i] = true
            Zarr[j] = true
            Xarr[i] = true
            Xarr[j] = true
        elseif direction == :Z
            Zarr[i] = true
            Zarr[j] = true
        end
        operators[k] = PauliOperator(0x00, Xarr, Zarr)
    end
end

function _HC_WilsonPlaquette_operator(i::Int64, L::Int64)
    i = i #X
    j = _HC_zneighbour(i, L) #Y
    k = _HC_xneighbour(j, L) #Z
    l = _HC_yneighbour(k, L) #X
    m = _HC_zneighbour(l, L) #Y
    n = _HC_xneighbour(m, L) #Z
    Xarr = falses(2*L^2)
    Zarr = falses(2*L^2)

    Xarr[i] = true

    Xarr[j] = true
    Zarr[j] = true

    Zarr[k] = true

    Xarr[l] = true

    Xarr[m] = true
    Zarr[m] = true

    Zarr[n] = true

    return PauliOperator(0x00, Xarr, Zarr)
end

function _HC_WilsonPlaquette_operators(L::Int64)
    operators = Vector{PauliOperator}(undef, L^2-1)
    for i in eachindex(operators)
        operators[i] = _HC_WilsonPlaquette_operator(2*i, L)
    end
    return operators
end

function _HC_WilsonLoops(L::Int64)
    site1 = 1
    site2 = 1
    Xarr1 = falses(2*L^2)
    Zarr1 = falses(2*L^2)
    Xarr2 = falses(2*L^2)
    Zarr2 = falses(2*L^2)
    for i in 1:2*L
        Xarr1[site1] = true #XXXX...
        Zarr2[site2] = true #ZZZZ...
        if isodd(i)
            site1 = _HC_zneighbour(site1, L)
            site2 = _HC_xneighbour(site2, L)
        else
            site1 = _HC_yneighbour(site1, L)
            site2 = _HC_yneighbour(site2, L)
        end
    end
    return [PauliOperator(0x00, Xarr1, Zarr1) PauliOperator(0x00, Xarr2, Zarr2)]
end

function _HC_interaction_operators(L::Int64)
end