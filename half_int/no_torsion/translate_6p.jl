function translate_6p_vec_to_ind(index::Int)

    temp = index

    i = 0
    U1 = 0
    U2 = 0
    U3 = 0
    F = 0

    while temp > 0

        i += 1

        temp -= len_i_3U_F[i]

    end

    temp += len_i_3U_F[i]

    while temp > 0

        U1 += 1

        temp -= len_3U_F[i,i,U1]

    end

    temp += len_3U_F[i,i,U1]

    (a1,b1) = SU[i,i,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_2U_F[a1,b1,U2]

    end

    temp += len_2U_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_U_F[a2,b2,U3]

    end

    temp += len_U_F[a2,b2,U3]

    F = temp

    return (i,U1,U2,U3,F)

end


function translate_6p_vec_to_ind_naive(index::Int)

    count = 0

    i0 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    F0 = 0

    for i in 1:k+1

        for U1 in 1:RangeU[i,i]

            (i2,a1,b1) = SU[i,i,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i3,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i4,a3,b3) = SU[a2,b2,U3,:]

                    for F in 1:RangeF[a3,b3]

                        (i5,i6) = SF[a3,b3,F,:]

                        count += 1

                        if count == index

                            (i0,U10,U20,U30,F0) =
                                (i,U1,U2,U3,F)

                        end

                    end

                end

            end

        end

    end

    return (i0,U10,U20,U30,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_6p_ind_to_vec(i,U1,U2,U3,F)

    temp = 0

    if i > 1

        temp += sum(len_i_3U_F[1:i-1])

    end

    if U1 > 1

        temp += sum(len_3U_F[i,i,1:U1-1])

    end

    (a1,b1) = SU[i,i,U1,2:3]

    if U2 > 1

        temp += sum(len_2U_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_U_F[a2,b2,1:U3-1])

    end

    temp += F

    return temp

end


#Test cases

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,U11,U21,U31,F1) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,U11,U21,U31,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,U11,U21,U31,F1) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,U11,U21,U31,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,U11,U21,U31,F1) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,U11,U21,U31,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,U11,U21,U31,F1) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,U11,U21,U31,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_vec_to_ind(test1))
#println(@time translate_6p_vec_to_ind_naive(test1))
#(i1,U11,U21,U31,F1) = translate_6p_vec_to_ind(test1)
#println(@time translate_6p_ind_to_vec(i1,U11,U21,U31,F1) - test1)
