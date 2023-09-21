function _chain_XX(size::Int) ::Vector{PauliOperator}
    XX = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Xarr[i] = true
        Xarr[mod1(i+1, size)] = true
        XX[i] = PauliOperator(0x00, Xarr, falses(size))
    end
    return unique(XX)
end

function _chain_YY(size::Int) ::Vector{PauliOperator}
    YY = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Xarr[i] = true
        Xarr[mod1(i+1, size)] = true
        Zarr = falses(size)
        Zarr[i] = true
        Zarr[mod1(i+1, size)] = true
        YY[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(YY)
end

function _chain_ZZ(size::Int) ::Vector{PauliOperator}
    ZZ = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Zarr = falses(size)
        Zarr[i] = true
        Zarr[mod1(i+1, size)] = true
        ZZ[i] = PauliOperator(0x00, falses(size), Zarr)
    end
    return unique(ZZ)
end

function _chain_XY(size::Int) ::Vector{PauliOperator}
    XY = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Xarr[i] = true
        Xarr[mod1(i+1, size)] = true
        Zarr[mod1(i+1, size)] = true
        XY[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(XY)
end

function _chain_XZ(size::Int) ::Vector{PauliOperator}
    XZ = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Xarr[i] = true
        Zarr[mod1(i+1, size)] = true
        XZ[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(XZ)
end

function _chain_YX(size::Int) ::Vector{PauliOperator}
    YX = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Xarr[i] = true
        Zarr[i] = true
        Xarr[mod1(i+1, size)] = true
        YX[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(YX)
end

function _chain_YZ(size::Int) ::Vector{PauliOperator}
    YZ = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Zarr[i] = true
        Zarr[mod1(i+1, size)] = true
        Xarr[i] = true
        YZ[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(YZ)
end

function _chain_ZX(size::Int) ::Vector{PauliOperator}
    ZX = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Zarr[i] = true
        Xarr[mod1(i+1, size)] = true
        ZX[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(ZX)
end

function _chain_ZY(size::Int) ::Vector{PauliOperator}
    ZY = Vector{PauliOperator}(undef, size)
    for i in 1:size
        Xarr = falses(size)
        Zarr = falses(size)
        Zarr[i] = true
        Xarr[mod1(i+1, size)] = true
        Zarr[mod1(i+1, size)] = true
        ZY[i] = PauliOperator(0x00, Xarr, Zarr)
    end
    return unique(ZY)
end

function _chain_blue(size::Int) #ZZ
    blue = Vector{PauliOperator}(undef, Int(size/3))
    for i in 1:size
        if mod1(i, 3) == 1
            Xarr = falses(size)
            Zarr = falses(size)
            Zarr[i] = true
            Zarr[mod1(i+1, size)] = true
            blue[div(i,3)+1] = PauliOperator(0x00, Xarr, Zarr)
        end
    end
    return blue
end

function _chain_red(size::Int) #XX
    red = Vector{PauliOperator}(undef, Int(size/3))
    for i in 1:size
        if mod1(i, 3) == 2
            Xarr = falses(size)
            Zarr = falses(size)
            Xarr[i] = true
            Xarr[mod1(i+1, size)] = true
            red[div(i,3)+1] = PauliOperator(0x00, Xarr, Zarr)
        end
    end
    return red
end

function _chain_green(size::Int) #YY
    green = Vector{PauliOperator}(undef, Int(size/3))
    for i in 1:size
        if mod1(i, 3) == 2
            Xarr = falses(size)
            Zarr = falses(size)
            Xarr[i] = true
            Xarr[mod1(i+1, size)] = true
            Zarr[i] = true
            Zarr[mod1(i+1, size)] = true
            green[div(i,3)+1] = PauliOperator(0x00, Xarr, Zarr)
        end
    end
    return green
end

function random_operator(trajectory::PPChainTrajectory, operators) ::PauliOperator
    probability = rand()

    if probability < trajectory.params[1]
        return operators[1, rand(1:trajectory.size)]
    elseif probability < (trajectory.params[1] + trajectory.params[2])
        return operators[2, rand(1:trajectory.size)]
    else
        return operators[3, rand(1:trajectory.size)]
    end
end
    
function random_operator(trajectory::PQChainTrajectory, operators) ::PauliOperator
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]

    probability1 = rand()
    probability2 = rand()

    if probability1 < px
        if probability2 < px
            return operators[1, rand(1:trajectory.size)]
        elseif probability2 < (px + py)
            return operators[2, rand(1:trajectory.size)]
        else
            return operators[3, rand(1:trajectory.size)]
        end
    elseif probability1 < (px + py)
        if probability2 < px
            return operators[4, rand(1:trajectory.size)]
        elseif probability2 < (px + py)
            return operators[5, rand(1:trajectory.size)]
        else
            return operators[6, rand(1:trajectory.size)]
        end
    else
        if probability2 < px
            return operators[7, rand(1:trajectory.size)]
        elseif probability2 < (px + py)
            return operators[8, rand(1:trajectory.size)]
        else
            return operators[9, rand(1:trajectory.size)]
        end
    end
end

function random_operator(trajectory::KekuleChainTrajectory, operators) ::PauliOperator
    probability = rand()

    if probability < trajectory.params[1]
        return operators[1, rand(1:Int(trajectory.size/3))]
    elseif probability < (trajectory.params[1] + trajectory.params[2])
        return operators[2, rand(1:Int(trajectory.size/3))]
    else
        return operators[3, rand(1:Int(trajectory.size/3))]
    end
end

function get_operators(trajectory::PPChainTrajectory)
    XX = _chain_XX(trajectory.size)
    YY = _chain_YY(trajectory.size)
    ZZ = _chain_ZZ(trajectory.size)
    matrix = Matrix{PauliOperator}(undef, 3, trajectory.size)
    for i in 1:trajectory.size
        matrix[1, i] = XX[i]
        matrix[2, i] = YY[i]
        matrix[3, i] = ZZ[i]
    end
    return matrix
end

function get_operators(trajectory::PQChainTrajectory)
    XX = _chain_XX(trajectory.size)
    XY = _chain_XY(trajectory.size)
    XZ = _chain_XZ(trajectory.size)
    YX = _chain_YX(trajectory.size)
    YY = _chain_YY(trajectory.size)
    YZ = _chain_YZ(trajectory.size)
    ZX = _chain_ZX(trajectory.size)
    ZY = _chain_ZY(trajectory.size)
    ZZ = _chain_ZZ(trajectory.size)
    matrix = Matrix{PauliOperator}(undef, 9, trajectory.size)
    for i in 1:trajectory.size
        matrix[1, i] = XX[i]
        matrix[2, i] = XY[i]
        matrix[3, i] = XZ[i]
        matrix[4, i] = YX[i]
        matrix[5, i] = YY[i]
        matrix[6, i] = YZ[i]
        matrix[7, i] = ZX[i]
        matrix[8, i] = ZY[i]
        matrix[9, i] = ZZ[i]
    end
    return matrix
end

function get_operators(trajectory::KekuleChainTrajectory)
    red = _chain_red(trajectory.size)
    green = _chain_green(trajectory.size)
    blue = _chain_blue(trajectory.size)
    matrix = Matrix{PauliOperator}(undef, 3, Int(trajectory.size/3))
    for i in 1:trajectory.size
        matrix[1, i] = red[i]
        matrix[2, i] = green[i]
        matrix[3, i] = blue[i]
    end
    return matrix
end