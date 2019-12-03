#Define the necessary lengths for the 10p states:
#First we can start from the len_4U_E

function length_5U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_4U[a1,b1,:])

        end

    end

    return len

end

len_5U = length_5U();


function length_6U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_5U[a1,b1,:])

        end

    end

    return len

end

len_6U = length_6U();


function length_7U()

    len7U = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len7U[a,b,U] = sum(len_6U[a1,b1,:])

        end

    end

    return len7U

end

len_7U = length_7U();


function length_8U()

    len8U = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len8U[a,b,U] = sum(len_7U[a1,b1,:])

        end

    end

    return len8U

end

len_8U = length_8U();


function length_ij_8U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_8U[i,j,:])

    end

    return len

end

len_ij_8U = length_ij_8U();

#println(len_s_K_8U_E)

configs_10p = sum(len_ij_8U)

#println(configs_10p * 16 / 1024 / 1024 / 1024)

#Define more indices for other 10p state basis

function length_X_2U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for X in 1:RangeX[a,b,i,j]

            (a1,b1) = SX[a,b,i,j,X,:]

            len[a,b,i,j,X] = sum(len_2U[a1,b1,:])

        end

    end

    return len

end

len_X_2U = length_X_2U();

function length_U_X_2U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,3:4]

            len[a,b,i,j,U] = sum(len_X_2U[a,b,a1,b1,:])

        end

    end

    return len

end

len_U_X_2U = length_U_X_2U();


function length_ij_U_X_2U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        len[a,b,i,j] = sum(len_U_X_2U[a,b,i,j,:])

    end

    return len

end

len_ij_U_X_2U = length_ij_U_X_2U();


function length_U_ij_U_X_2U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_U_X_2U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_2U = length_U_ij_U_X_2U();


function length_2U_ij_U_X_2U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_U_ij_U_X_2U[a1,b1,:])

        end

    end

    return len

end

len_2U_ij_U_X_2U = length_2U_ij_U_X_2U();


function length_3U_ij_U_X_2U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_2U_ij_U_X_2U[a1,b1,:])

        end

    end

    return len

end

len_3U_ij_U_X_2U = length_3U_ij_U_X_2U();


function length_4U_ij_U_X_2U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_3U_ij_U_X_2U[a1,b1,:])

        end

    end

    return len

end

len_4U_ij_U_X_2U = length_4U_ij_U_X_2U();


function length_ij_4U_ij_U_X_2U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_4U_ij_U_X_2U[i,j,:])

    end

    return len

end

len_ij_4U_ij_U_X_2U = length_ij_4U_ij_U_X_2U();

#println(len_s_K_4U_s_K_U_X_2U_E)

#println(sum(len_ij_4U_ij_U_X_2U[:])) # same as configs_10p

#println(configs_10p)

function length_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for X in 1:RangeX[a,b,a1,b1]

            (a2,b2) = SX[a,b,a1,b1,X,:]

            len[a,b,a1,b1,X] = RangeU[a2,b2]

        end

    end

    return len

end

len_X_U = length_X_U();


function length_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a,b,a1,b1,U] = sum(len_X_U[a,b,a2,b2,:])

        end

    end

    return len

end

len_U_X_U = length_U_X_U();


function length_2U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,3:4]

            len[a,b,a1,b1,U] = sum(len_U_X_U[a,b,a2,b2,:])

        end

    end

    return len

end

len_2U_X_U = length_2U_X_U();


function length_ij_2U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        len[a,b,i,j] = sum(len_2U_X_U[a,b,i,j,:])

    end

    return len

end

len_ij_2U_X_U = length_ij_2U_X_U();



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

len_U_ij_2U_X_U = length_U_ij_2U_X_U();


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

len_2U_ij_2U_X_U = length_2U_ij_2U_X_U();


function length_3U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_2U_ij_2U_X_U[a1,b1,:])

        end

    end

    return len

end

len_3U_ij_2U_X_U = length_3U_ij_2U_X_U();


function length_4U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_3U_ij_2U_X_U[a1,b1,:])

        end

    end

    return len

end

len_4U_ij_2U_X_U = length_4U_ij_2U_X_U();


function length_ij_4U_ij_2U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_4U_ij_2U_X_U[i,j,:])

    end

    return len

end

len_ij_4U_ij_2U_X_U = length_ij_4U_ij_2U_X_U(); #same number of configs as 10p.

#println(sum(len_ij_4U_ij_2U_X_U))

#println(configs_10p)

# Lengths for the 4th 10p transformation

# can use the ones for U_X_U_E!

function length_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i in 1:k+1, j in 1:k+1

            len[a,b,i,j] = sum(len_U_X_U[a,b,i,j,:])

        end

    end

    return len

end

len_ij_U_X_U = length_ij_U_X_U();


function length_U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_U = length_U_ij_U_X_U();


function length_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_U_ij_U_X_U[a1,b1,:])

        end

    end

    return len

end

len_2U_ij_U_X_U = length_2U_ij_U_X_U();



function length_3U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_2U_ij_U_X_U[a1,b1,:])

        end

    end

    return len

end

len_3U_ij_U_X_U = length_3U_ij_U_X_U();



function length_4U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_3U_ij_U_X_U[a1,b1,:])

        end

    end

    return len

end

len_4U_ij_U_X_U = length_4U_ij_U_X_U();



function length_5U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_4U_ij_U_X_U[a1,b1,:])

        end

    end

    return len

end

len_5U_ij_U_X_U = length_5U_ij_U_X_U();




function length_ij_5U_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_5U_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_5U_ij_U_X_U = length_ij_5U_ij_U_X_U();

#println(sum(len_ij_5U_ij_U_X_U))

#println(configs_10p)

#Define the lengths necessary to define the 10 p state before svd

# We can use everything up to len_X_U_E again.

function length_X_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x2 in 1:k+1, y2 in 1:k+1

        for X in 1:RangeX[a,b,x2,y2]

            (a2,b2) = SX[a,b,x2,y2,X,:]

            len[a,b,x2,y2,X] = sum(len_2U_ij_U_X_U[a2,b2,:])

        end

    end

    return len

end

len_X_2U_ij_U_X_U = length_X_2U_ij_U_X_U();


function length_U_X_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (x2,y2) = SU[a1,b1,U,3:4]

            len[a,b,a1,b1,U] = sum(len_X_2U_ij_U_X_U[a,b,x2,y2,:])

        end

    end

    return len

end

len_U_X_2U_ij_U_X_U = length_U_X_2U_ij_U_X_U();


function length_ij_U_X_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i2 in 1:k+1, j2 in 1:k+1

        len[a,b,i2,j2] = sum(len_U_X_2U_ij_U_X_U[a,b,i2,j2,:])

    end

    return len

end

len_ij_U_X_2U_ij_U_X_U = length_ij_U_X_2U_ij_U_X_U();


function length_U_ij_U_X_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_ij_U_X_2U_ij_U_X_U[a1,b1,:,:])

        end

    end

    return len

end

len_U_ij_U_X_2U_ij_U_X_U = length_U_ij_U_X_2U_ij_U_X_U();


function length_ij_U_ij_U_X_2U_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_U_ij_U_X_2U_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_U_ij_U_X_2U_ij_U_X_U = length_ij_U_ij_U_X_2U_ij_U_X_U();

#configs_10p_bSVD = sum(len_s_K_U_s_K_U_E_X_2U_s_K_U_E_X_U_E[:])
#println(sum(len_ij_U_ij_U_X_2U_ij_U_X_U[:,:]))

#println(configs_10p)
