#Ingredients necessary to vectorize the 6p states.

#Here I compute and store the lengths of indices for the 6p states.
#This is necessary to quickly calculate the basis states from the vector index.


function length_2U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = RangeU[a1,b1]

        end

    end

    return len

end

len_2U = length_2U();


function length_3U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_2U[a1,b1,:])

        end

    end

    return len

end

len_3U = length_3U();


function length_4U()

    len = zeros(Int,k+1,k+1,maxU)

    for a in 1:k+1, b in 1:k+1

        for U in 1:RangeU[a,b]

            (a1,b1) = SU[a,b,U,3:4]

            len[a,b,U] = sum(len_3U[a1,b1,:])

        end

    end

    return len

end

len_4U = length_4U();


function length_ij_4U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_4U[i,j,:])

    end

    return len

end

len_ij_4U = length_ij_4U();

#Total length of 6p state:

configs_6p = sum(len_ij_4U)

#println(configs_6p)
