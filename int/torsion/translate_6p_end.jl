function translate_6p_end_vec_to_ind(index::Int)
    #s5,K5,U1,s6,K6,U2,X,U3,E

    temp = index

    P5 = 0
    U1 = 0
    P6 = 0
    U2 = 0
    X = 0
    U3 = 0

    while temp > 0

        P5 += 1

        (i5,j5) = SP[P5,:]

        temp -= len_ij_U_ij_U_X_U[i5,j5]

    end

    (i5,j5) = SP[P5,:]

    temp += len_ij_U_ij_U_X_U[i5,j5]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_U[i5,j5,U1]

    end

    temp += len_U_ij_U_X_U[i5,j5,U1]

    (a1,b1) = SU[i5,j5,U1,3:4]

    while temp > 0

        P6 += 1

        (i6,j6) = SP[P6,:]

        temp -= len_ij_U_X_U[a1,b1,i6,j6]

    end

    (i6,j6) = SP[P6,:]

    temp += len_ij_U_X_U[a1,b1,i6,j6]

    while temp > 0

        U2 += 1

        temp -= len_U_X_U[a1,b1,i6,j6,U2]

    end

    temp += len_U_X_U[a1,b1,i6,j6,U2]

    (a2,b2) = SU[i6,j6,U2,3:4]

    while temp > 0

        X += 1

        temp -= len_X_U[a1,b1,a2,b2,X]

    end

    temp += len_X_U[a1,b1,a2,b2,X]

    (a3,b3) = SX[a1,b1,a2,b2,X,:]

    U3 = temp

    return (i5,j5,U1,i6,j6,U2,X,U3)

end


function translate_6p_end_vec_to_ind_naive(index::Int)

    count = 0

    i50 = 0
    j50 = 0
    U10 = 0
    i60 = 0
    j60 = 0
    U20 = 0
    X0 = 0
    U30 = 0

    for i5 in 1:k+1, j5 in 1:k+1

        for U1 in 1:RangeU[i5,j5]

            (i1,j1,a1,b1) = SU[i5,j5,U1,:]

            for i6 in 1:k+1, j6 in 1:k+1

                for U2 in 1:RangeU[i6,j6]

                    (i3,j3,a4,b4) = SU[i6,j6,U2,:]

                    for X in 1:RangeX[a1,b1,a4,b4]

                        (a3,b3) = SX[a1,b1,a4,b4,X,:]

                        for U3 in 1:RangeU[a3,b3]

                            (i2,j2,i4,j4) = SU[a3,b3,U3,:]

                            count += 1

                            if count == index

                                    (i50,j50,U10,i60,j60,U20,X0,U30) =
                                        (i5,j5,U1,i6,j6,U2,X,U3)

                            end

                        end

                    end

                end

            end

        end

    end

    return (i50,j50,U10,i60,j60,U20,X0,U30)

end

#Define the translation from set of indices to the vectorized index.

function translate_6p_end_ind_to_vec(i5::Int,j5::Int,U1::Int,
    i6::Int,j6::Int,U2::Int,X::Int,
    U3::Int)

    temp = 0

    P5 = Pf[i5,j5]

    if P5 > 1

        for P5ind in 1:P5-1

            (i5ind,j5ind) = SP[P5ind,:]

            temp += len_ij_U_ij_U_X_U[i5ind,j5ind]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_U[i5,j5,1:U1-1])

    end

    (a1,b1) = SU[i5,j5,U1,3:4]

    P6 = Pf[i6,j6]

    if P6 > 1

        for P6ind in 1:P6-1

            (i6ind,j6ind) = SP[P6ind,:]

            temp += len_ij_U_X_U[a1,b1,i6ind,j6ind]

        end

    end

    if U2 > 1

        temp += sum(len_U_X_U[a1,b1,i6,j6,1:U2-1])

    end

    (a2,b2) = SU[i6,j6,U2,3:4]

    if X > 1

        temp += sum(len_X_U[a1,b1,a2,b2,1:X-1])

    end

    (a3,b3) = SX[a1,b1,a2,b2,X,:]

    temp += U3

    return temp

end


#Test cases

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,j51,U11,i61,j61,U21,X1,U31) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,j51,U11,i61,j61,U21,X1,U31) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,j51,U11,i61,j61,U21,X1,U31) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,j51,U11,i61,j61,U21,X1,U31) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,j51,U11,i61,j61,U21,X1,U31) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,j51,U11,i61,j61,U21,X1,U31) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,j51,U11,i61,j61,U21,X1,U31) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,j51,U11,i61,j61,U21,X1,U31) - test1)

#test1 = rand(1:configs_6p)
#println(@time translate_6p_end_vec_to_ind(test1))
#println(@time translate_6p_end_vec_to_ind_naive(test1))
#(i51,j51,U11,i61,j61,U21,X1,U31) = translate_6p_end_vec_to_ind(test1)
#println(@time translate_6p_end_ind_to_vec(i51,j51,U11,i61,j61,U21,X1,U31) - test1)
