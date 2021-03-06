function length_X_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1,maxX)

    for a1 in 1:k+1, b1 in 1:k+1, x1 in 1:k+1, y1 in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        for X in 1:RangeX[a1,b1,x1,y1]

            (a2,b2) = SX[a1,b1,x1,y1,X,:]

            len[a1,b1,x1,y1,x3,y3,X] = sum(len_X_U[a2,b2,x3,y3,:])

        end

    end

    return len

end

len_X_X_U_svd = length_X_X_U_svd();


function length_U_X_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, i1 in 1:k+1, j1 in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        for U in 1:RangeU[i1,j1]

            (x1,y1) = SU[i1,j1,U,3:4]

            len[a1,b1,i1,j1,x3,y3,U] = sum(len_X_X_U_svd[a1,b1,x1,y1,x3,y3,:])

        end

    end

    return len

end

len_U_X_X_U_svd = length_U_X_X_U_svd();


function length_ij_U_X_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1)

    for a1 in 1:k+1, b1 in 1:k+1, i1 in 1:k+1, j1 in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        len[a1,b1,i1,j1,x3,y3] = sum(len_U_X_X_U_svd[a1,b1,i1,j1,x3,y3,:])

    end

    return len

end

len_ij_U_X_X_U_svd = length_ij_U_X_X_U_svd();


function length_U_ij_U_X_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a1,b1,x3,y3,U] = sum(len_ij_U_X_X_U_svd[a2,b2,:,:,x3,y3])

        end

    end

    return len

end

len_U_ij_U_X_X_U_svd = length_U_ij_U_X_X_U_svd();


function length_ij_U_ij_U_X_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for x2 in 1:k+1, y2 in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        len[x2,y2,x3,y3] = sum(len_U_ij_U_X_X_U_svd[x2,y2,x3,y3,:])

    end

    return len

end

len_ij_U_ij_U_X_X_U_svd = length_ij_U_ij_U_X_X_U_svd();

function length_ind2_8p_svd()

    len = zeros(Int,k+1,k+1)

    for x3 in 1:k+1, y3 in 1:k+1

        len[x3,y3] = sum(len_ij_U_ij_U_X_X_U_svd[:,:,x3,y3])

    end

    return len

end

configs_8p_ind2_svd = length_ind2_8p_svd()

function length_ind1_8p_svd()

    len = zeros(Int,k+1,k+1)

    for x3 in 1:k+1, y3 in 1:k+1

        len[x3,y3] = RangeU[x3,y3]

    end

    return len

end

configs_8p_ind1_svd = length_ind1_8p_svd()

function compute_8p_max_ind1()

    temp = 0

    for x3 in 1:k+1, y3 in 1:k+1

        temp1 = configs_8p_ind1_svd[x3,y3]

        if temp1 > temp

            temp = temp1

        end

    end

    return temp

end

maxind1_8p = compute_8p_max_ind1();

function compute_8p_max_ind2()

    temp = 0

    for x3 in 1:k+1, y3 in 1:k+1

        temp1 = configs_8p_ind2_svd[x3,y3]

        if temp1 > temp

            temp = temp1

        end

    end

    return temp

end

maxind2_8p = compute_8p_max_ind2();
