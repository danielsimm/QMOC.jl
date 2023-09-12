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

function random_operator(trajectory::PPChainTrajectory) ::PauliOperator
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]

    probability = rand()

    if probability < px
        return rand(trajectory.XX_projectors)
    elseif probability < (px + py)
        return rand(trajectory.YY_projectors)
    else
        return rand(trajectory.ZZ_projectors)
    end
end
    
function random_operator(trajectory::PQChainTrajectory) ::PauliOperator
    px = trajectory.params[1]
    py = trajectory.params[2]
    pz = trajectory.params[3]

    probability1 = rand()
    probability2 = rand()

    if probability1 < px
        if probability2 < px
            return rand(trajectory.XX_projectors)
        elseif probability2 < (px + py)
            return rand(trajectory.XY_projectors)
        else
            return rand(trajectory.XZ_projectors)
        end
    elseif probability1 < (px + py)
        if probability2 < px
            return rand(trajectory.YX_projectors)
        elseif probability2 < (px + py)
            return rand(trajectory.YY_projectors)
        else
            return rand(trajectory.YZ_projectors)
        end
    else
        if probability2 < px
            return rand(trajectory.ZX_projectors)
        elseif probability2 < (px + py)
            return rand(trajectory.ZY_projectors)
        else
            return rand(trajectory.ZZ_projectors)
        end
    end
end