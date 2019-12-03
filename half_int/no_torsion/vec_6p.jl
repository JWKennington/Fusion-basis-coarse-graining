#Ingredients necessary to vectorize the 6p states.

#Here I compute and store the lengths of indices for the 6p states.
#This is necessary to quickly calculate the basis states from the vector index.

function length_U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = RangeF[a1,b1]

        end

    end

    return len

end

len_U_F = length_U_F();


function length_2U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,2:3]

            len[a,b,U] = sum(len_U_F[a1,b1,:])

        end

    end

    return len

end

len_2U_F = length_2U_F();


function length_3U_F()

    len = zeros(Int,k+1,k+1,maxU)

    for i in 1:k+1, j in 1:k+1

        for U in 1:RangeU[i,j]

            (a1,b1) = SU[i,j,U,2:3]

            len[i,j,U] = sum(len_2U_F[a1,b1,:])

        end

    end

    return len

end

len_3U_F = length_3U_F();



function length_i_3U_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_3U_F[i,i,:])

    end

    return len

end

len_i_3U_F = length_i_3U_F();

#println(len_i_3U_F)

#Total length of 6p state:

configs_6p = sum(len_i_3U_F[:])

#println(configs_6p)
