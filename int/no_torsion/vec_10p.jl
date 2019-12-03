#Define the necessary lengths for the 10p states:
#First we can start from the len_3U_F


function length_4U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_3U_F[a1,b1,:])

        end

    end

    return len

end

len_4U_F = length_4U_F();


function length_5U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_4U_F[a1,b1,:])

        end

    end

    return len

end

len_5U_F = length_5U_F();


function length_6U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_5U_F[a1,b1,:])

        end

    end

    return len

end

len_6U_F = length_6U_F();


function length_7U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_6U_F[a1,b1,:])

        end

    end

    return len

end

len_7U_F = length_7U_F();


function length_i_7U_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_7U_F[i,i,:])

    end

    return len

end

len_i_7U_F = length_i_7U_F();

configs_10p = sum(len_i_7U_F)

#println(configs_10p)


#Define more indices for other 10p state basis

function length_X_U_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for X in 1:RangeX[a,b,i,j]

            (a1,b1) = SX[a,b,i,j,X,:]

            len[a,b,i,j,X] = sum(len_U_F[a1,b1,:])

        end

    end

    return len

end

len_X_U_F = length_X_U_F();


function length_U_X_U_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[a,b,i,j,U] = sum(len_X_U_F[a,b,a1,b1,:])

        end

    end

    return len

end

len_U_X_U_F = length_U_X_U_F();

function length_i_U_X_U_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1

        len[a,b,i] = sum(len_U_X_U_F[a,b,i,i,:])

    end

    return len

end

len_i_U_X_U_F = length_i_U_X_U_F();


function length_U_i_U_X_U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_i_U_X_U_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_U_F = length_U_i_U_X_U_F();



function length_2U_i_U_X_U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_U_i_U_X_U_F[a1,b1,:])

        end

    end

    return len

end

len_2U_i_U_X_U_F = length_2U_i_U_X_U_F();


function length_3U_i_U_X_U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_2U_i_U_X_U_F[a1,b1,:])

        end

    end

    return len

end

len_3U_i_U_X_U_F = length_3U_i_U_X_U_F();


function length_4U_i_U_X_U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_3U_i_U_X_U_F[a1,b1,:])

        end

    end

    return len

end

len_4U_i_U_X_U_F = length_4U_i_U_X_U_F();



function length_i_4U_i_U_X_U_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_4U_i_U_X_U_F[i,i,:])

    end

    return len

end

len_i_4U_i_U_X_U_F = length_i_4U_i_U_X_U_F(); #total number of configs same as configs_10p

#println(len_i_4U_i_U_X_U_F)

#println(sum(len_i_4U_i_U_X_U_F))

#println(configs_10p)

function length_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for X in 1:RangeX[a,b,a1,b1]

            (a2,b2) = SX[a,b,a1,b1,X,:]

            len[a,b,a1,b1,X] = RangeF[a2,b2]

        end

    end

    return len

end

len_X_F = length_X_F();


function length_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, a1 in 1:k+1, b1 in 1:k+1

        for U in 1:RangeU[a1,b1]

            (a2,b2) = SU[a1,b1,U,2:3]

            len[a,b,a1,b1,U] = sum(len_X_F[a,b,a2,b2,:])

        end

    end

    return len

end

len_U_X_F = length_U_X_F();



function length_2U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a2,b2) = SU[i,j,U,2:3]

            len[a,b,i,j,U] = sum(len_U_X_F[a,b,a2,b2,:])

        end

    end

    return len

end

len_2U_X_F = length_2U_X_F();



function length_i_2U_X_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1

        len[a,b,i] = sum(len_2U_X_F[a,b,i,i,:])

    end

    return len

end

len_i_2U_X_F = length_i_2U_X_F();



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

len_U_i_2U_X_F = length_U_i_2U_X_F();



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

len_2U_i_2U_X_F = length_2U_i_2U_X_F();



function length_3U_i_2U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_2U_i_2U_X_F[a1,b1,:])

        end

    end

    return len

end

len_3U_i_2U_X_F = length_3U_i_2U_X_F();



function length_4U_i_2U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_3U_i_2U_X_F[a1,b1,:])

        end

    end

    return len

end

len_4U_i_2U_X_F = length_4U_i_2U_X_F();



function length_i_4U_i_2U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_4U_i_2U_X_F[i,i,:])

    end

    return len

end

len_i_4U_i_2U_X_F = length_i_4U_i_2U_X_F(); #same number of configs as 10p.

#println(len_i_4U_i_2U_X_F)

#println(sum(len_i_4U_i_2U_X_F))

#println(configs_10p)


# Lengths for the 4th 10p transformation



function length_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        for i in 1:k+1

            len[a,b,i] = sum(len_U_X_F[a,b,i,i,:])

        end

    end

    return len

end

len_i_U_X_F = length_i_U_X_F();



function length_U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_F = length_U_i_U_X_F();


function length_2U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_U_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_2U_i_U_X_F = length_2U_i_U_X_F();



function length_3U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_2U_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_3U_i_U_X_F = length_3U_i_U_X_F();



function length_4U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_3U_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_4U_i_U_X_F = length_4U_i_U_X_F();



function length_5U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_4U_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_5U_i_U_X_F = length_5U_i_U_X_F();




function length_i_5U_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_5U_i_U_X_F[i,i,:])

    end

    return len

end

len_i_5U_i_U_X_F = length_i_5U_i_U_X_F(); #same as configs_10p


#println(len_i_5U_i_U_X_F)

#println(sum(len_i_5U_i_U_X_F))

#println(configs_10p)


#Define the lengths necessary to define the 10 p state before svd

# We can use everything up to len_X_F again.

#function length_E_X_F()

#    len = zeros(Int,k+1,k+1,k+1,k+1,maxE)

#    for a in 1:k+1, b in 1:k+1, x4 in 1:k+1, y4 in 1:k+1

#        for E in 1:RangeE[x4,y4]

#            len[a,b,x4,y4,E] = sum(len_X_F[a,b,x4,y4,:])

#        end

#    end

#    return len

#end

#len_E_X_F = length_E_X_F();


#function length_U_E_X_F()

#    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

#    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

#        for U in 1:RangeU[i,j]

#            (x4,y4) = SU[i,j,U,2:3]

#            len[a,b,i,j,U] = sum(len_E_X_F[a,b,x4,y4,:])

#        end

#    end

#    return len

#end

#len_U_E_X_F = length_U_E_X_F();



#function length_i_U_E_X_F()

#    len = zeros(Int,k+1,k+1,k+1)

#    for a in 1:k+1, b in 1:k+1, i4 in 1:k+1

#        len[a,b,i4] = sum(len_U_E_X_F[a,b,i4,i4,:])

#    end

#    return len

#end

#len_i_U_E_X_F = length_i_U_E_X_F();



#function length_U_i_U_E_X_F()

#    len = zeros(Int,k+1,k+1,maxU)

#    for a in 1:k+1, b in 1:k+1

#        for U in 1:RangeU[a,b]

#            (a1,b1) = SU[a,b,U,2:3]

#            len[a,b,U] = sum(len_i_U_E_X_F[a1,b1,:])

#        end

#    end

#    return len

#end

#len_U_i_U_E_X_F = length_U_i_U_E_X_F();



#function length_2U_i_U_E_X_F()

#    len = zeros(Int,k+1,k+1,maxU)

#    for a in 1:k+1, b in 1:k+1

#        for U in 1:RangeU[a,b]

#            (a1,b1) = SU[a,b,U,2:3]

#            len[a,b,U] = sum(len_U_i_U_E_X_F[a1,b1,:])

#        end

#    end

#    return len

#end

#len_2U_i_U_E_X_F = length_2U_i_U_E_X_F();

#We have all the necesssary length already available.
# The last one: len_2U_i_U_X_F

function length_X_2U_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxX)

    for a in 1:k+1, b in 1:k+1, x2 in 1:k+1, y2 in 1:k+1

        for X in 1:RangeX[a,b,x2,y2]

            (a2,b2) = SX[a,b,x2,y2,X,:]

            len[a,b,x2,y2,X] = sum(len_2U_i_U_X_F[a2,b2,:])

        end

    end

    return len

end

len_X_2U_i_U_X_F = length_X_2U_i_U_X_F();

#The following we do not need any more. E not necessary

#function length_E_X_2U_i_U_E_X_F()

#    len = zeros(Int,k+1,k+1,k+1,k+1,maxE)

#    for a in 1:k+1, b in 1:k+1, x2 in 1:k+1, y2 in 1:k+1

#        for E in 1:RangeE[x2,y2]

#            len[a,b,x2,y2,E] = sum(len_X_2U_i_U_E_X_F[a,b,x2,y2,:])

#        end

#    end

#    return len

#end

#len_E_X_2U_i_U_E_X_F = length_E_X_2U_i_U_E_X_F();



function length_U_X_2U_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1, i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (x2,y2) = SU[i,j,U,2:3]

            len[a,b,i,j,U] = sum(len_X_2U_i_U_X_F[a,b,x2,y2,:])

        end

    end

    return len

end

len_U_X_2U_i_U_X_F = length_U_X_2U_i_U_X_F();




function length_i_U_X_2U_i_U_X_F()

    len = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1, i2 in 1:k+1

        len[a,b,i2] = sum(len_U_X_2U_i_U_X_F[a,b,i2,i2,:])

    end

    return len

end

len_i_U_X_2U_i_U_X_F = length_i_U_X_2U_i_U_X_F();



function length_U_i_U_X_2U_i_U_X_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_i_U_X_2U_i_U_X_F[a1,b1,:])

        end

    end

    return len

end

len_U_i_U_X_2U_i_U_X_F = length_U_i_U_X_2U_i_U_X_F();




function length_i_U_i_U_X_2U_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_U_i_U_X_2U_i_U_X_F[i,i,:])

    end

    return len

end

len_i_U_i_U_X_2U_i_U_X_F = length_i_U_i_U_X_2U_i_U_X_F();


#println(configs_10p*16/1024/1024/1024)


#println(sum(len_i_U_i_U_X_2U_i_U_X_F)) #same as configs_10p

#println(configs_10p)

#println(configs_10p_bSVD)
