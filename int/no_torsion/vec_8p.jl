#Here I define all the lengths necessary for the 8p states

#lengths for the 8p state after SVD.


function length_i_5U_F()

    len = zeros(Int,k+1)

    for i1 in 1:k+1

        len[i1] = sum(len_5U_F[i1,i1,:])

    end

    return len

end

len_i_5U_F = length_i_5U_F();

configs_8p = sum(len_i_5U_F[:])


#Lengths for trafo1:

function length_X_2U_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x in 1:k+1, y in 1:k+1

        for X1 in 1:RangeX[a,b,x,y]

            (a2,b2) = SX[a,b,x,y,X1,:]

            len[a,b,x,y,X1] = sum(len_2U_F[a2,b2,:])

        end

    end

    return len

end

len_X_2U_F = length_X_2U_F();


function length_U_X_2U_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (x,y) = SU[i,j,U,2:3]

            len[a,b,i,j,U] = sum(len_X_2U_F[a,b,x,y,:])

        end

    end

    return len

end

len_U_X_2U_F = length_U_X_2U_F();



function length_i_U_X_2U_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1

        len[a,b,i] = sum(len_U_X_2U_F[a,b,i,i,:])

    end

    return len

end

len_i_U_X_2U_F = length_i_U_X_2U_F();



function length_U_i_U_X_2U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i2 in 1:k+1, j2 in 1:k+1

        for U in 1:RangeU[i2,j2]

            (a1,b1) = SU[i2,j2,U,2:3]

            len[i2,j2,U] = sum(len_i_U_X_2U_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_2U_F = length_U_i_U_X_2U_F();




function length_i_U_i_U_X_2U_F()

    len = zeros(Int,k+1)

    for u2 in 1:k+1

        len[u2] = sum(len_U_i_U_X_2U_F[u2,u2,:])

    end

    return len

end

len_i_U_i_U_X_2U_F = length_i_U_i_U_X_2U_F();

#Total amount of configs identical to configs_8p.

#Lengths necessary for trafo_2

function length_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x in 1:k+1, y in 1:k+1

        for X in 1:RangeX[a,b,x,y]

            (a1,b1) = SX[a,b,x,y,X,:]

            len[a,b,x,y,X] = sum(len_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_X_i_U_X_F = length_X_i_U_X_F();



function length_U_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (x,y) = SU[i,j,U,2:3]

            len[a,b,i,j,U] = sum(len_X_i_U_X_F[a,b,x,y,:])

        end

    end

    return len

end

len_U_X_i_U_X_F = length_U_X_i_U_X_F();




function length_i_U_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i in 1:k+1

            len[a,b,i] = sum(len_U_X_i_U_X_F[a,b,i,i,:])

        end

    end

    return len

end

len_i_U_X_i_U_X_F = length_i_U_X_i_U_X_F();



function length_U_i_U_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_i_U_X_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_i_U_X_F = length_U_i_U_X_i_U_X_F();



function length_i_U_i_U_X_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_U_i_U_X_i_U_X_F[i,i,:])

    end

    return len

end

len_i_U_i_U_X_i_U_X_F = length_i_U_i_U_X_i_U_X_F();

# Total number of configs agrees with configs_8p

# Lengths for trafo_3


function length_i_3U_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_3U_i_U_X_F[i,i,:])

    end

    return len

end

len_i_3U_i_U_X_F = length_i_3U_i_U_X_F()

# Again number of total configs agrees with configs_8p.

# Lengths for trafo_4

function length_U_i_2U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_i_2U_X_F[a1,b1,:])

        end

    end

    return len

end


len_U_i_2U_X_F = length_U_i_2U_X_F()


function length_2U_i_2U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_U_i_2U_X_F[a1,b1,:])

        end

    end

    return len

end

len_2U_i_2U_X_F = length_2U_i_2U_X_F()




function length_i_2U_i_2U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_2U_i_2U_X_F[i,i,:])

    end

    return len

end

len_i_2U_i_2U_X_F = length_i_2U_i_2U_X_F()


# Number of total configs matches configs_8p

# Define the necessary lengths for trafo_6 (trafo_5 uses the same indices as trafo_4)
# trafo_6 is the same as trafo_3

#Now define the lengths for the state before the svd


function length_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x3 in 1:k+1, y3 in 1:k+1

        for X in 1:RangeX[a,b,x3,y3]

            (a1,b1) = SX[a,b,x3,y3,X,:]

            len[a,b,x3,y3,X] = sum(len_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_X_i_U_X_F = length_X_i_U_X_F();




function length_U_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i1 in 1:k+1, j1 in 1:k+1

        for U in 1:RangeU[i1,j1]

            (x1,y1) = SU[i1,j1,U,2:3]

            len[a,b,i1,j1,U] = sum(len_X_i_U_X_F[a,b,x1,y1,:])

        end

    end

    return len

end

len_U_X_i_U_X_F = length_U_X_i_U_X_F();



function length_i_E_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i1 in 1:k+1

            len[a,b,i1] = sum(len_U_X_i_U_X_F[a,b,i1,i1,:])

        end

    end

    return len

end

len_i_U_X_i_U_X_F = length_i_U_X_i_U_X_F();



function length_U_i_U_X_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_i_U_X_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_i_U_X_F = length_U_i_U_X_i_U_X_F();




function length_i_U_i_U_X_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_U_i_U_X_i_U_X_F[i,i,:])

    end

    return len

end

len_i_U_i_U_X_i_U_X_F = length_i_U_i_U_X_i_U_X_F();

#println(sum(len_i_U_i_U_X_i_U_X_F[:])) #same as configs_8p

#println(configs_8p)
