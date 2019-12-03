#Define the necessary length to define the matrix
#that will be split into two via an SVD

#First the smaller part that only has two indices, U, E. Already defined!

#Larger part, need new indices. Defined up to X_U_E. But, here is a difference.
#We have to explicitly keep the entries of X from the top, since these are the
#ones kept fixed during an SVD.


function length_U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a1,b1,x4,y4,U] = sum(len_X_U[a2,b2,x4,y4,:])

        end

    end

    return len

end

len_U_X_U_svd = length_U_X_U_svd();


function length_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a1,b1,x4,y4,U] = sum(len_U_X_U_svd[a2,b2,x4,y4,:])

        end

    end

    return len

end

len_2U_X_U_svd = length_2U_X_U_svd();


function length_X_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1,maxX)

    for a1 in 1:k+1, b1 in 1:k+1, x2 in 1:k+1, y2 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        for X in 1:RangeX[a1,b1,x2,y2]

            (a2,b2) = SX[a1,b1,x2,y2,X,:]

            len[a1,b1,x2,y2,x4,y4,X] = sum(len_2U_X_U_svd[a2,b2,x4,y4,:])

        end

    end

    return len

end

len_X_2U_X_U_svd = length_X_2U_X_U_svd();


function length_U_X_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, i2 in 1:k+1, j2 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        for U in 1:RangeU[i2,j2]

            (x2,y2) = SU[i2,j2,U,3:4]

            len[a1,b1,i2,j2,x4,y4,U] = sum(len_X_2U_X_U_svd[a1,b1,x2,y2,x4,y4,:])

        end

    end

    return len

end

len_U_X_2U_X_U_svd = length_U_X_2U_X_U_svd();


function length_ij_U_X_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1)

    for a1 in 1:k+1, b1 in 1:k+1, i2 in 1:k+1, j2 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        len[a1,b1,i2,j2,x4,y4] = sum(len_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4,:])

    end

    return len

end

len_ij_U_X_2U_X_U_svd = length_ij_U_X_2U_X_U_svd();


function length_U_ij_U_X_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a1 in 1:k+1, b1 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a1,b1,x4,y4,U] = sum(len_ij_U_X_2U_X_U_svd[a2,b2,:,:,x4,y4])

        end

    end

    return len

end

len_U_ij_U_X_2U_X_U_svd = length_U_ij_U_X_2U_X_U_svd();


function length_ij_U_ij_U_X_2U_X_U_svd()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for i1 in 1:k+1, j1 in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

        len[i1,j1,x4,y4] = sum(len_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4,:])

    end

    return len

end

len_ij_U_ij_U_X_2U_X_U_svd = length_ij_U_ij_U_X_2U_X_U_svd();

function length_ind2_10p_svd()

    len = zeros(Int,k+1,k+1)

    for x4 in 1:k+1, y4 in 1:k+1

        len[x4,y4] = sum(len_ij_U_ij_U_X_2U_X_U_svd[:,:,x4,y4])

    end

    return len

end

configs_10p_ind2_svd = length_ind2_10p_svd()

function length_ind1_10p_svd()

    len = zeros(Int,k+1,k+1)

    for x4 in 1:k+1, y4 in 1:k+1

        len[x4,y4] = RangeU[x4,y4]

    end

    return len

end

configs_10p_ind1_svd = length_ind1_10p_svd()

function compute_max_ind1()

    temp = 0

    for x4 in 1:k+1, y4 in 1:k+1

        temp1 = configs_10p_ind1_svd[x4,y4]

        if temp1 > temp

            temp = temp1

        end

    end

    return temp

end

maxind1 = compute_max_ind1();

function compute_max_ind2()

    temp = 0

    for x4 in 1:k+1, y4 in 1:k+1

        temp1 = configs_10p_ind2_svd[x4,y4]

        if temp1 > temp

            temp = temp1

        end

    end

    return temp

end

maxind2 = compute_max_ind2();

#println(configs_10p_ind1_svd)
#println(configs_10p_ind2_svd)
