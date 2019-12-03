function translate_8p_SVD_ind2_vec_to_ind(index::Int,x3::Int,y3::Int)

    #u2,K2,U1,s1,K1,U2,E1p,X1,X2,U4,E

    i2 = 0
    U1 = 0
    i1 = 0
    U2 = 0
    X1 = 0
    X2 = 0
    F = 0

    temp = index

    while temp > 0

        i2 += 1

        temp -= len_i_U_i_U_X_X_F_svd[i2,x3,y3]

    end

    temp += len_i_U_i_U_X_X_F_svd[i2,x3,y3]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_X_F_svd[i2,i2,x3,y3,U1]

    end

    temp += len_U_i_U_X_X_F_svd[i2,i2,x3,y3,U1]

    (a1,b1) = SU[i2,i2,U1,2:3]

    while temp > 0

        i1 += 1

        temp -= len_i_U_X_X_F_svd[a1,b1,i1,x3,y3]

    end

    temp += len_i_U_X_X_F_svd[a1,b1,i1,x3,y3]

    while temp > 0

        U2 += 1

        temp -= len_U_X_X_F_svd[a1,b1,i1,i1,x3,y3,U2]

    end

    temp += len_U_X_X_F_svd[a1,b1,i1,i1,x3,y3,U2]

    (x1,y1) = SU[i1,i1,U2,2:3]

    while temp > 0

        X1 += 1

        temp -= len_X_X_F_svd[a1,b1,x1,y1,x3,y3,X1]

    end

    temp += len_X_X_F_svd[a1,b1,x1,y1,x3,y3,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        X2 += 1

        temp -= len_X_F[a2,b2,x3,y3,X2]

    end

    temp += len_X_F[a2,b2,x3,y3,X2]

    F = temp

    return (i2,U1,i1,U2,X1,X2,F)

end



function translate_8p_SVD_ind2_vec_to_ind_naive(index::Int,x3::Int,y3::Int)

    i20 = 0
    U10 = 0
    i10 = 0
    U20 = 0
    X10 = 0
    X20 = 0
    F0 = 0

    count = 0

    for i2 in 1:k+1

        for U1 in 1:RangeU[i2,i2]

            (i5,M1,N1) = SU[i2,i2,U1,:]

            for i1 in 1:k+1

                for U2 in 1:RangeU[i1,i1]

                    (k1,x1,y1) = SU[i1,i1,U2,:]

                    for X1 in 1:RangeX[M1,N1,x1,y1]

                        (M5,N5) = SX[M1,N1,x1,y1,X1,:]

                        for X2 in 1:RangeX[M5,N5,x3,y3]

                            (M3,N3) = SX[M5,N5,x3,y3,X2,:]

                            for F in 1:RangeF[M3,N3]

                                (k6,x4) = SF[M3,N3,F,:]

                                count += 1

                                if count == index

                                    (i20,U10,i10,U20,X10,X20,F0) =
                                        (i2,U1,i1,U2,X1,X2,F)

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,U10,i10,U20,X10,X20,F0)

end



function translate_8p_SVD_ind2_ind_to_vec(i2::Int,U1::Int,
    i1::Int,U2::Int,X1::Int,X2::Int,F::Int,
    x3::Int,y3::Int)

    temp = 0

    if i2 > 1

        temp += sum(len_i_U_i_U_X_X_F_svd[1:i2-1,x3,y3])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_X_F_svd[i2,i2,x3,y3,1:U1-1])

    end

    (a1,b1) = SU[i2,i2,U1,2:3]

    if i1 > 1

        temp += sum(len_i_U_X_X_F_svd[a1,b1,1:i1-1,x3,y3])

    end

    if U2 > 1

        temp += sum(len_U_X_X_F_svd[a1,b1,i1,i1,x3,y3,1:U2-1])

    end

    (x1,y1) = SU[i1,i1,U2,2:3]

    if X1 > 1

        temp += sum(len_X_X_F_svd[a1,b1,x1,y1,x3,y3,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if X2 > 1

        temp += sum(len_X_F[a2,b2,x3,y3,1:X2-1])

    end

    temp += F

    return temp

end

function testing_8p(x3::Int,y3::Int)

    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i2,U1,i1,U2,X1,X2,F) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i2,U1,i1,U2,X1,X2,F,x3,y3) - test1)


    test2 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test2,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test2,x3,y3))
    (i2,U1,i1,U2,X1,X2,F) = translate_8p_SVD_ind2_vec_to_ind(test2,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i2,U1,i1,U2,X1,X2,F,x3,y3) - test2)


    test3 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test3,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test3,x3,y3))
    (i2,U1,i1,U2,X1,X2,F) = translate_8p_SVD_ind2_vec_to_ind(test3,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i2,U1,i1,U2,X1,X2,F,x3,y3) - test3)


    test4 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test4,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test4,x3,y3))
    (i2,U1,i1,U2,X1,X2,F) = translate_8p_SVD_ind2_vec_to_ind(test4,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i2,U1,i1,U2,X1,X2,F,x3,y3) - test4)


    test5 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test5,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test5,x3,y3))
    (i2,U1,i1,U2,X1,X2,F) = translate_8p_SVD_ind2_vec_to_ind(test5,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i2,U1,i1,U2,X1,X2,F,x3,y3) - test5)

end

#testing_8p(1,1)
#testing_8p(1,2)
#testing_8p(2,1)
#testing_8p(2,2)

#No need here, since the index is just F.

#function translate_8p_SVD_ind1_vec_to_ind(index::Int,x3::Int,y3::Int)

#    Z3 = 0
#    E3 = 0

#    temp = index

#    while temp > 0

#        Z3 += 1

#        temp -= len_U_E[x3,y3,Z3]

#    end

#    temp += len_U_E[x3,y3,Z3]

#    E3 = temp

#    return (Z3,E3)

#end

#function translate_8p_SVD_ind1_vec_to_ind_naive(index::Int,x3::Int,y3::Int)

#    Z30 = 0
#    E30 = 0

#    count = 0

#    for Z3 in 1:RangeU[x3,y3]

#        (a1,b1) = SU[x3,y3,Z3,3:4]

#        for E3 in 1:RangeE[a1,b1]

#            count += 1

#            if count == index

#                (Z30,E30) = (Z3,E3)

#            end

#        end

#    end

#    return (Z30,E30)

#end


#function translate_8p_SVD_ind1_ind_to_vec(Z3::Int,E3::Int,x3::Int,y3::Int)

#    temp = 0

#    if Z3 > 1

#        temp += sum(len_U_E[x3,y3,1:Z3-1])

#    end

#    temp += E3

#    return temp

#end

#function testing_8p_2(x3::Int,y3::Int)

#    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
#    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
#    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
#    (Z31,E31) = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
#    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,E31,x3,y3) - test1)


#    test2 = rand(1:configs_8p_ind1_svd[x3,y3])
#    println(@time translate_8p_SVD_ind1_vec_to_ind(test2,x3,y3))
#    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test2,x3,y3))
#    (Z32,E32) = translate_8p_SVD_ind1_vec_to_ind(test2,x3,y3)
#    println(@time translate_8p_SVD_ind1_ind_to_vec(Z32,E32,x3,y3) - test2)


#    test3 = rand(1:configs_8p_ind1_svd[x3,y3])
#    println(@time translate_8p_SVD_ind1_vec_to_ind(test3,x3,y3))
#    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test3,x3,y3))
#    (Z33,E33) = translate_8p_SVD_ind1_vec_to_ind(test3,x3,y3)
#    println(@time translate_8p_SVD_ind1_ind_to_vec(Z33,E33,x3,y3) - test3)


#    test4 = rand(1:configs_8p_ind1_svd[x3,y3])
#    println(@time translate_8p_SVD_ind1_vec_to_ind(test4,x3,y3))
#    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test4,x3,y3))
#    (Z34,E34) = translate_8p_SVD_ind1_vec_to_ind(test4,x3,y3)
#    println(@time translate_8p_SVD_ind1_ind_to_vec(Z34,E34,x3,y3) - test4)


#    test5 = rand(1:configs_8p_ind1_svd[x3,y3])
#    println(@time translate_8p_SVD_ind1_vec_to_ind(test5,x3,y3))
#    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test5,x3,y3))
#    (Z35,E35) = translate_8p_SVD_ind1_vec_to_ind(test5,x3,y3)
#    println(@time translate_8p_SVD_ind1_ind_to_vec(Z35,E35,x3,y3) - test5)

#end

#testing_2(1,2)
