function _HC_all_partitions(L)
    partitions_per_direction = div(L, 4)
    partition_size = div(2*L^2, 4)
    xrows = _HC_xrows(L)
    yrows = _HC_yrows(L)
    zrows = _HC_zrows(L)
    partitions = []

    # z partitions
        for i in 1:partitions_per_direction
            A_rows = mod1.(collect(1:div(L, 4)) .+ (i-1), L)
            B_rows = mod1.(collect(div(L, 4)+1:div(L, 2)) .+ (i-1), L)
            C_rows = mod1.(collect(div(L, 2)+1:div(3L, 4)) .+ (i-1), L)
            D_rows = mod1.(collect(div(3L, 4)+1:L) .+ (i-1), L)
            push!(partitions, [reduce(vcat, zrows[A_rows]), reduce(vcat, zrows[B_rows]), reduce(vcat, zrows[C_rows]), reduce(vcat, zrows[D_rows])])
        end

    # y partitions
        for i in 1:partitions_per_direction
            A_rows = mod1.(collect(1:div(L, 4)) .+ (i-1), L)
            B_rows = mod1.(collect(div(L, 4)+1:div(L, 2)) .+ (i-1), L)
            C_rows = mod1.(collect(div(L, 2)+1:div(3L, 4)) .+ (i-1), L)
            D_rows = mod1.(collect(div(3L, 4)+1:L) .+ (i-1), L)
            push!(partitions, [reduce(vcat, yrows[A_rows]), reduce(vcat, yrows[B_rows]), reduce(vcat, yrows[C_rows]), reduce(vcat, yrows[D_rows])])
        end

    # x partitions
        for i in 1:partitions_per_direction
            A_rows = mod1.(collect(1:div(L, 4)) .+ (i-1), L)
            B_rows = mod1.(collect(div(L, 4)+1:div(L, 2)) .+ (i-1), L)
            C_rows = mod1.(collect(div(L, 2)+1:div(3L, 4)) .+ (i-1), L)
            D_rows = mod1.(collect(div(3L, 4)+1:L) .+ (i-1), L)
            push!(partitions, [reduce(vcat, xrows[A_rows]), reduce(vcat, xrows[B_rows]), reduce(vcat, xrows[C_rows]), reduce(vcat, xrows[D_rows])])
        end
        return partitions
end

function _HC_zrows(L)
    return [collect((2*L*(i-1) + 1):2*L*i) for i in 1:L]
end
function _HC_yrows(L)
    rows = [zeros(Int, 2*L) for _ in 1:L]
    for i in 1:L
        rows[i][1] = 2*(i-1) + 1
        for j in 2:2*L
            if iseven(j)
                rows[i][j] = rows[i][j-1] + 1
            else
                rows[i][j] = rows[i][j-1] + 2*L - 1
            end
        end
    end 
    return rows
end
function _HC_xrows(L)
    rows = [zeros(Int, 2*L) for _ in 1:L]
    for i in 1:L
        rows[i][1] = 2*L*(i-1) + 1
        for j in 2:2*L
            if isodd(j)
                rows[i][j] = QMOC._HC_yneighbour(rows[i][j-1], L)
            else
                rows[i][j] = QMOC._HC_zneighbour(rows[i][j-1], L)
            end
        end
    end 
    return rows
end