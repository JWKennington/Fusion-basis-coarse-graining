function translate_6p_vec_to_ind(index::Int)

    temp = index

    P = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0

    while temp > 0

        P += 1

        (i1,j1) = SP[P,:]

        temp -= len_ij_4U[i1,j1]

    end

    (i1,j1) = SP[P,:]

    temp += len_ij_4U[i1,j1]

    while temp > 0

        U1 += 1

        temp -= len_4U[i1,j1,U1]

    end

    temp += len_4U[i1,j1,U1]

    (a1,b1) = SU[i1,j1,U1,3:4]

    while temp > 0

        U2 += 1

        temp -= len_3U[a1,b1,U2]

    end

    temp += len_3U[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,3:4]

    while temp > 0

        U3 += 1

        temp -= len_2U[a2,b2,U3]

    end

    temp += len_2U[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,3:4]

    #println(temp)

    U4 = temp

    return (i1,j1,U1,U2,U3,U4)

end


function translate_6p_vec_to_ind_naive(index::Int)

    count = 0

    i10 = 0
    j10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0

    for i1 in 1:k+1, j1 in 1:k+1

        for U1 in 1:RangeU[i1,j1]

            (i2,j2,a1,b1) = SU[i1,j1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i3,j3,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i4,j4,a3,b3) = SU[a2,b2,U3,:]

                    for U4 in 1:RangeU[a3,b3]

                        (i5,j5,i6,j6) = SU[a3,b3,U4,:]

                        count += 1

                        if count == index

                            (i10,j10,U10,U20,U30,U40) =
                                (i1,j1,U1,U2,U3,U4)

                        end

                    end

                end

            end

        end

    end

    return (i10,j10,U10,U20,U30,U40)

end

#Define the translation from set of indices to the vectorized index.

function translate_6p_ind_to_vec(i::Int,j::Int,U1::Int,U2::Int,
    U3::Int,U4::Int)

    temp = 0

    P = Pf[i,j]

    if P > 1

        for Pind in 1:P-1

            (iind,jind) = SP[Pind,:]

            temp += len_ij_4U[iind,jind]

        end

    end

    if U1 > 1

        temp += sum(len_4U[i,j,1:U1-1])

    end

    (a1,b1) = SU[i,j,U1,3:4]

    if U2 > 1

        temp += sum(len_3U[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,3:4]

    if U3 > 1

        temp += sum(len_2U[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,3:4]

    temp += U4

    return temp

end


#Test cases

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,j1,U11,U21,U31,U41) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,j1,U11,U21,U31,U41) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,j1,U11,U21,U31,U41) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,j1,U11,U21,U31,U41) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,j1,U11,U21,U31,U41) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,j1,U11,U21,U31,U41) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,j1,U11,U21,U31,U41) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,j1,U11,U21,U31,U41) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,j1,U11,U21,U31,U41) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,j1,U11,U21,U31,U41) - test1)
