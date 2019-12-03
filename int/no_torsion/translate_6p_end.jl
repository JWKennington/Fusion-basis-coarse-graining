function translate_6p_end_vec_to_ind(index::Int)

    temp = index

    i5 = 0
    U1 = 0
    i6 = 0
    U2 = 0
    X = 0
    F = 0

    while temp > 0

        i5 += 1

        temp -= len_i_U_i_U_X_F[i5]

    end

    temp += len_i_U_i_U_X_F[i5]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_F[i5,i5,U1]

    end

    temp += len_U_i_U_X_F[i5,i5,U1]

    (a1,b1) = SU[i5,i5,U1,2:3]

    while temp > 0

        i6 += 1

        temp -= len_i_U_X_F[a1,b1,i6]

    end

    temp += len_i_U_X_F[a1,b1,i6]

    while temp > 0

        U2 += 1

        temp -= len_U_X_F[a1,b1,i6,i6,U2]

    end

    temp += len_U_X_F[a1,b1,i6,i6,U2]

    (a2,b2) = SU[i6,i6,U2,2:3]

    while temp > 0

        X += 1

        temp -= len_X_F[a1,b1,a2,b2,X]

    end

    temp += len_X_F[a1,b1,a2,b2,X]

    F = temp

    return (i5,U1,i6,U2,X,F)

end


function translate_6p_end_vec_to_ind_naive(index::Int)

    count = 0

    i50 = 0
    U10 = 0
    i60 = 0
    U20 = 0
    X0 = 0
    F0 = 0

    for i5 in 1:k+1

        for U1 in 1:RangeU[i5,i5]

            (i1,a1,b1) = SU[i5,i5,U1,:]

            for i6 in 1:k+1

                for U2 in 1:RangeU[i6,i6]

                    (i3,a4,b4) = SU[i6,i6,U2,:]

                    for X in 1:RangeX[a1,b1,a4,b4]

                        (a3,b3) = SX[a1,b1,a4,b4,X,:]

                        for F in 1:RangeF[a3,b3]

                            (i2,i4) = SU[a3,b3,F,:]

                            count += 1

                            if count == index

                                (i50,U10,i60,U20,X0,F0) =
                                    (i5,U1,i6,U2,X,F)

                            end

                        end

                    end

                end

            end

        end

    end

    return (i50,U10,i60,U20,X0,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_6p_end_ind_to_vec(i5::Int,U1::Int,
    i6::Int,U2::Int,X::Int,F::Int)

    temp = 0

    if i5 > 1

        temp += sum(len_i_U_i_U_X_F[1:i5-1])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_F[i5,i5,1:U1-1])

    end

    (a1,b1) = SU[i5,i5,U1,2:3]

    if i6 > 1

        temp += sum(len_i_U_X_F[a1,b1,1:i6-1])

    end

    if U2 > 1

        temp += sum(len_U_X_F[a1,b1,i6,i6,1:U2-1])

    end

    (a2,b2) = SU[i6,i6,U2,2:3]

    if X > 1

        temp += sum(len_X_F[a1,b1,a2,b2,1:X-1])

    end

    (a3,b3) = SX[a1,b1,a2,b2,X,:]

    temp += F

    return temp

end


#Test cases

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,U11,i61,U21,X1,F1) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,U11,i61,U21,X1,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,U11,i61,U21,X1,F1) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,U11,i61,U21,X1,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,U11,i61,U21,X1,F1) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,U11,i61,U21,X1,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,U11,i61,U21,X1,F1) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,U11,i61,U21,X1,F1) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,U11,i61,U21,X1,F1) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,U11,i61,U21,X1,F1) - test1)
