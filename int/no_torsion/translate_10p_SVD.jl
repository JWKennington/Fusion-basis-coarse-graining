function translate_10p_SVD_ind2_vec_to_ind(index::Int,x4::Int,y4::Int)

    i1 = 0
    U1 = 0
    i2 = 0
    Z2 = 0
    X2 = 0
    U4 = 0
    U5 = 0
    X4 = 0
    F = 0

    temp = index

    while temp > 0

        i1 += 1

        temp -= len_i_U_i_U_X_2U_X_F_svd[i1,x4,y4]

    end

    temp += len_i_U_i_U_X_2U_X_F_svd[i1,x4,y4]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_2U_X_F_svd[i1,x4,y4,U1]

    end

    temp += len_U_i_U_X_2U_X_F_svd[i1,x4,y4,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        i2 += 1

        temp -= len_i_U_X_2U_X_F_svd[a1,b1,i2,x4,y4]

    end

    temp += len_i_U_X_2U_X_F_svd[a1,b1,i2,x4,y4]

    while temp > 0

        Z2 += 1

        temp -= len_U_X_2U_X_F_svd[a1,b1,i2,x4,y4,Z2]

    end

    temp += len_U_X_2U_X_F_svd[a1,b1,i2,x4,y4,Z2]

    (x2,y2) = SU[i2,i2,Z2,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_2U_X_F_svd[a1,b1,x2,y2,x4,y4,X2]

    end

    temp += len_X_2U_X_F_svd[a1,b1,x2,y2,x4,y4,X2]

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    while temp > 0

        U4 += 1

        temp -= len_2U_X_F_svd[a2,b2,x4,y4,U4]

    end

    temp += len_2U_X_F_svd[a2,b2,x4,y4,U4]

    (a3,b3) = SU[a2,b2,U4,2:3]

    while temp > 0

        U5 += 1

        temp -= len_U_X_F_svd[a3,b3,x4,y4,U5]

    end

    temp += len_U_X_F_svd[a3,b3,x4,y4,U5]

    (a4,b4) = SU[a3,b3,U5,2:3]

    while temp > 0

        X4 += 1

        temp -= len_X_F[a4,b4,x4,y4,X4]

    end

    temp += len_X_F[a4,b4,x4,y4,X4]

    F = temp

    return (i1,U1,i2,Z2,X2,U4,U5,X4,F)

end



function translate_10p_SVD_ind2_vec_to_ind_naive(index::Int,x4::Int,y4::Int)

    i10 = 0
    U10 = 0
    i20 = 0
    Z20 = 0
    X20 = 0
    U40 = 0
    U50 = 0
    X40 = 0
    F0 = 0

    count = 0

    for i1 in 1:k+1

        for U1 in 1:RangeU[i1,i1]

            (i3,a1,b1) = SU[i1,i1,U1,:]

            for i2 in 1:k+1

                for Z2 in 1:RangeU[i2,i2]

                    (k2,x2,y2) = SU[i2,i2,Z2,:]

                    for X2 in 1:RangeX[a1,b1,x2,y2]

                        (a3,b3) = SX[a1,b1,x2,y2,X2,:]

                        for U4 in 1:RangeU[a3,b3]

                            (i5,M2,N2) = SU[a3,b3,U4,:]

                            for U5 in 1:RangeU[M2,N2]

                                (k6,c2,d2) = SU[M2,N2,U5,:]

                                for X4 in 1:RangeX[c2,d2,x4,y4]

                                    (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                                    for F in 1:RangeF[c1,d1]

                                        count += 1

                                        if count == index

                                            (i10,U10,i20,Z20,X20,U40,U50,X40,F0) =
                                                (i1,U1,i2,Z2,X2,U4,U5,X4,F)

                                        end

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i10,U10,i20,Z20,X20,U40,U50,X40,F0)

end



function translate_10p_SVD_ind2_ind_to_vec(i1::Int,U1::Int,
    i2::Int,Z2::Int,X2::Int,
    U4::Int,U5::Int,X4::Int,F::Int,
    x4::Int,y4::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_U_i_U_X_2U_X_F_svd[1:i1-1,x4,y4])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_2U_X_F_svd[i1,x4,y4,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if i2 > 1

        temp += sum(len_i_U_X_2U_X_F_svd[a1,b1,1:i2-1,x4,y4])

    end

    if Z2 > 1

        temp += sum(len_U_X_2U_X_F_svd[a1,b1,i2,x4,y4,1:Z2-1])

    end

    (x2,y2) = SU[i2,i2,Z2,2:3]

    if X2 > 1

        temp += sum(len_X_2U_X_F_svd[a1,b1,x2,y2,x4,y4,1:X2-1])

    end

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    if U4 > 1

        temp += sum(len_2U_X_F_svd[a2,b2,x4,y4,1:U4-1])

    end

    (a3,b3) = SU[a2,b2,U4,2:3]

    if U5 > 1

        temp += sum(len_U_X_F_svd[a3,b3,x4,y4,1:U5-1])

    end

    (a4,b4) = SU[a3,b3,U5,2:3]

    if X4 > 1

        temp += sum(len_X_F[a4,b4,x4,y4,1:X4-1])

    end

    temp += F

    return temp

end

function testing(x4::Int,y4::Int)

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,U11,i21,Z21,X21,U41,U51,X41,F1) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,U11,i21,Z21,X21,U41,U51,X41,F1,x4,y4) - test1)

    println()

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,U11,i21,Z21,X21,U41,U51,X41,F1) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,U11,i21,Z21,X21,U41,U51,X41,F1,x4,y4) - test1)

    println()

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,U11,i21,Z21,X21,U41,U51,X41,F1) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,U11,i21,Z21,X21,U41,U51,X41,F1,x4,y4) - test1)

    println()

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,U11,i21,Z21,X21,U41,U51,X41,F1) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,U11,i21,Z21,X21,U41,U51,X41,F1,x4,y4) - test1)

    println()

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,U11,i21,Z21,X21,U41,U51,X41,F1) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,U11,i21,Z21,X21,U41,U51,X41,F1,x4,y4) - test1)

end

#testing(1,1)
#testing(1,2)
#testing(2,1)
#testing(2,2)

#No need, because we do not allow the "old" punctures to carry torsion, we only need index F.

#function translate_10p_SVD_ind1_vec_to_ind(index::Int,x4::Int,y4::Int)

    #Z4 = 0
    #E4 = 0

    #temp = index

    #while temp > 0

    #    Z4 += 1

    #    temp -= len_U_E[x4,y4,Z4]

    #end

    #temp += len_U_E[x4,y4,Z4]

    #E4 = temp

    #return (Z4,E4)

#end

#function translate_10p_SVD_ind1_vec_to_ind_naive(index::Int,x4::Int,y4::Int)

    #Z40 = 0
    #E40 = 0

    #count = 0

    #for Z4 in 1:RangeU[x4,y4]

    #    (a1,b1) = SU[x4,y4,Z4,3:4]

    #    for E4 in 1:RangeE[a1,b1]

    #        count += 1

    #        if count == index

    #            (Z40,E40) = (Z4,E4)

    #        end

    #    end

    #end

    #return (Z40,E40)

#end


#function translate_10p_SVD_ind1_ind_to_vec(Z4::Int,E4::Int,x4::Int,y4::Int)

    #temp = 0

    #if Z4 > 1

    #    temp += sum(len_U_E[x4,y4,1:Z4-1])

    #end

    #temp += E4

    #return temp

#end

#function testing_2(x4::Int,y4::Int)

    #test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    #println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    #println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    #(Z41,E41) = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    #println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,E41,x4,y4) - test1)


    #test2 = rand(1:configs_10p_ind1_svd[x4,y4])
    #println(@time translate_10p_SVD_ind1_vec_to_ind(test2,x4,y4))
    #println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test2,x4,y4))
    #(Z42,E42) = translate_10p_SVD_ind1_vec_to_ind(test2,x4,y4)
    #println(@time translate_10p_SVD_ind1_ind_to_vec(Z42,E42,x4,y4) - test2)


    #test3 = rand(1:configs_10p_ind1_svd[x4,y4])
    #println(@time translate_10p_SVD_ind1_vec_to_ind(test3,x4,y4))
    #println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test3,x4,y4))
    #(Z43,E43) = translate_10p_SVD_ind1_vec_to_ind(test3,x4,y4)
    #println(@time translate_10p_SVD_ind1_ind_to_vec(Z43,E43,x4,y4) - test3)


    #test4 = rand(1:configs_10p_ind1_svd[x4,y4])
    #println(@time translate_10p_SVD_ind1_vec_to_ind(test4,x4,y4))
    #println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test4,x4,y4))
    #(Z44,E44) = translate_10p_SVD_ind1_vec_to_ind(test4,x4,y4)
    #println(@time translate_10p_SVD_ind1_ind_to_vec(Z44,E44,x4,y4) - test4)


    #test5 = rand(1:configs_10p_ind1_svd[x4,y4])
    #println(@time translate_10p_SVD_ind1_vec_to_ind(test5,x4,y4))
    #println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test5,x4,y4))
    #(Z45,E45) = translate_10p_SVD_ind1_vec_to_ind(test5,x4,y4)
    #println(@time translate_10p_SVD_ind1_ind_to_vec(Z45,E45,x4,y4) - test5)

#end

#testing_2(1,2)
