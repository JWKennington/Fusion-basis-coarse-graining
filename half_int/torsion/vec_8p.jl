#Here I define all the lengths necessary for the 8p states

#lengths for the 8p state after SVD.

function length_ij_6U()

    len = zeros(Int,k+1,k+1)

    for i1 in 1:k+1, j1 in 1:k+1

        len[i1,j1] = sum(len_6U[i1,j1,:])

    end

    return len

end

len_ij_6U = length_ij_6U();

configs_8p = sum(len_ij_6U[:,:])

#println(configs_8p)


#Lengths for trafo1:

function length_X_3U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x in 1:k+1, y in 1:k+1

        for X1 in 1:RangeX[a,b,x,y]

            (a2,b2) = SX[a,b,x,y,X1,:]

            len[a,b,x,y,X1] = sum(len_3U[a2,b2,:])

        end

    end

    return len

end

len_X_3U = length_X_3U();


function length_U_X_3U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (x,y) = SU[i,j,U,3:4]

            len[a,b,i,j,U] = sum(len_X_3U[a,b,x,y,:])

        end

    end

    return len

end

len_U_X_3U = length_U_X_3U();



function length_ij_U_X_3U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        len[a,b,i,j] = sum(len_U_X_3U[a,b,i,j,:])

    end

    return len

end

len_ij_U_X_3U = length_ij_U_X_3U();


function length_U_ij_U_X_3U()

    len = zeros(Int,k+1,k+1,maxU)

    for i2 in 1:k+1, j2 in 1:k+1

        for U in 1:RangeU[i2,j2]

            (a1,b1) = SU[i2,j2,U,3:4]

            len[i2,j2,U] = sum(len_ij_U_X_3U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_3U = length_U_ij_U_X_3U();


function length_ij_U_ij_U_X_3U()

    len = zeros(Int,k+1,k+1)

    for x2 in 1:k+1, y2 in 1:k+1

        len[x2,y2] = sum(len_U_ij_U_X_3U[x2,y2,:])

    end

    return len

end

len_ij_U_ij_U_X_3U = length_ij_U_ij_U_X_3U();

#println(sum(len_ij_U_ij_U_X_3U))

#Total amount of configs identical to configs_8p.

#Lengths necessary for trafo_2

function length_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x in 1:k+1, y in 1:k+1

        for X in 1:RangeX[a,b,x,y]

            (a1,b1) = SX[a,b,x,y,X,:]

            len[a,b,x,y,X] = sum(len_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_X_ij_U_X_U = length_X_ij_U_X_U();


function length_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (x,y) = SU[i,j,U,3:4]

            len[a,b,i,j,U] = sum(len_X_ij_U_X_U[a,b,x,y,:])

        end

    end

    return len

end

len_U_X_ij_U_X_U = length_U_X_ij_U_X_U();


function length_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i in 1:k+1, j in 1:k+1

            len[a,b,i,j] = sum(len_U_X_ij_U_X_U[a,b,i,j,:])

        end

    end

    return len

end

len_ij_U_X_ij_U_X_U = length_ij_U_X_ij_U_X_U();


function length_U_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_U_X_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_ij_U_X_U = length_U_ij_U_X_ij_U_X_U();


function length_ij_U_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_U_ij_U_X_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_U_ij_U_X_ij_U_X_U = length_ij_U_ij_U_X_ij_U_X_U();

#println(sum(len_ij_U_ij_U_X_ij_U_X_U))

# Total number of configs agrees with configs_8p

# Lengths for trafo_3

function length_ij_3U_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_3U_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_3U_ij_U_X_U = length_ij_3U_ij_U_X_U()

#println(sum(len_ij_3U_ij_U_X_U))

# Again number of total configs agrees with configs_8p.

# Lengths for trafo_4

function length_U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_2U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_2U_X_U = length_U_ij_2U_X_U()


function length_2U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_U_ij_2U_X_U[a1,b1,:])

        end

    end

    return len

end

len_2U_ij_2U_X_U = length_2U_ij_2U_X_U()


function length_ij_2U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_2U_ij_2U_X_U[i,j,:])

    end

    return len

end

len_ij_2U_ij_2U_X_U = length_ij_2U_ij_2U_X_U()

#println(sum(len_ij_2U_ij_2U_X_U))

# Number of total configs matches configs_8p

# Define the necessary lengths for trafo_6 (trafo_5 uses the same indices as trafo_4)
# trafo_6 is the same as trafo_3

#Now define the lengths for the state before the svd


function length_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        for X in 1:RangeX[a,b,x3,y3]

            (a1,b1) = SX[a,b,x3,y3,X,:]

            len[a,b,x3,y3,X] = sum(len_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_X_ij_U_X_U = length_X_ij_U_X_U();


function length_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i1 in 1:k+1, j1 in 1:k+1

        for U in 1:RangeU[i1,j1]

            (x1,y1) = SU[i1,j1,U,3:4]

            len[a,b,i1,j1,U] = sum(len_X_ij_U_X_U[a,b,x1,y1,:])

        end

    end

    return len

end

len_U_X_ij_U_X_U = length_U_X_ij_U_X_U();


function length_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i1 in 1:k+1, j1 in 1:k+1

            len[a,b,i1,j1] = sum(len_U_X_ij_U_X_U[a,b,i1,j1,:])

        end

    end

    return len

end

len_ij_U_X_ij_U_X_U = length_ij_U_X_ij_U_X_U();


function length_U_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_U_X_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_ij_U_X_U = length_U_ij_U_X_ij_U_X_U();


function length_ij_U_ij_U_X_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_U_ij_U_X_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_U_ij_U_X_ij_U_X_U = length_ij_U_ij_U_X_ij_U_X_U();

#println(sum(len_ij_U_ij_U_X_ij_U_X_U))

#configs_8p_bSVD = sum(len_s_K_U_s_K_U_E_X_s_K_U_E_X_U_E[:])
