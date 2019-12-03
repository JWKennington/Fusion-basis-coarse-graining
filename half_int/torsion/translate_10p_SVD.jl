function translate_10p_SVD_ind2_vec_to_ind(index::Int,x4::Int,y4::Int)

    P = 0
    U1 = 0
    P2 = 0
    Z2 = 0
    X2 = 0
    U4 = 0
    U5 = 0
    X4 = 0
    U8 = 0

    temp = index

    while temp > 0

        P += 1

        (i1,j1) = SP[P,:]

        temp -= len_ij_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4]

    end

    (i1,j1) = SP[P,:]

    temp += len_ij_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4,U1]

    end

    temp += len_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4,U1]

    (a1,b1) = SU[i1,j1,U1,3:4]

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4]

    while temp > 0

        Z2 += 1

        temp -= len_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4,Z2]

    end

    temp += len_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4,Z2]

    (x2,y2) = SU[i2,j2,Z2,3:4]

    while temp > 0

        X2 += 1

        temp -= len_X_2U_X_U_svd[a1,b1,x2,y2,x4,y4,X2]

    end

    temp += len_X_2U_X_U_svd[a1,b1,x2,y2,x4,y4,X2]

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    while temp > 0

        U4 += 1

        temp -= len_2U_X_U_svd[a2,b2,x4,y4,U4]

    end

    temp += len_2U_X_U_svd[a2,b2,x4,y4,U4]

    (a3,b3) = SU[a2,b2,U4,3:4]

    while temp > 0

        U5 += 1

        temp -= len_U_X_U_svd[a3,b3,x4,y4,U5]

    end

    temp += len_U_X_U_svd[a3,b3,x4,y4,U5]

    (a4,b4) = SU[a3,b3,U5,3:4]

    while temp > 0

        X4 += 1

        temp -= len_X_U[a4,b4,x4,y4,X4]

    end

    temp += len_X_U[a4,b4,x4,y4,X4]

    (a5,b5) = SX[a4,b4,x4,y4,X4,:]

    U8 = temp

    return (i1,j1,U1,i2,j2,Z2,X2,U4,U5,X4,U8)

end



function translate_10p_SVD_ind2_vec_to_ind_naive(index::Int,x4::Int,y4::Int)

    i10 = 0
    j10 = 0
    U10 = 0
    i20 = 0
    j20 = 0
    Z20 = 0
    X20 = 0
    U40 = 0
    U50 = 0
    X40 = 0
    U80 = 0

    count = 0

    for i1 in 1:k+1, j1 in 1:k+1

        for U1 in 1:RangeU[i1,j1]

            (i3,j3,a1,b1) = SU[i1,j1,U1,:]

            for i2 in 1:k+1, j2 in 1:k+1

                for Z2 in 1:RangeU[i2,j2]

                    (k2,l2,x2,y2) = SU[i2,j2,Z2,:]

                    for X2 in 1:RangeX[a1,b1,x2,y2]

                        (a3,b3) = SX[a1,b1,x2,y2,X2,:]

                        for U4 in 1:RangeU[a3,b3]

                            (i5,j5,M2,N2) = SU[a3,b3,U4,:]

                            for U5 in 1:RangeU[M2,N2]

                                (k6,l6,c2,d2) = SU[M2,N2,U5,:]

                                for X4 in 1:RangeX[c2,d2,x4,y4]

                                    (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                                    for U8 in 1:RangeU[c1,d1]

                                        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

                                        count += 1

                                        if count == index

                                                (i10,j10,U10,i20,j20,Z20,X20,U40,U50,X40,U80) =
                                                    (i1,j1,U1,i2,j2,Z2,X2,U4,U5,X4,U8)

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

    return (i10,j10,U10,i20,j20,Z20,X20,U40,U50,X40,U80)

end



function translate_10p_SVD_ind2_ind_to_vec(i1::Int,j1::Int,U1::Int,
    i2::Int,j2::Int,Z2::Int,X2::Int,
    U4::Int,U5::Int,X4::Int,U8::Int,
    x4::Int,y4::Int)

    temp = 0

    P = Pf[i1,j1]

    if P > 1

        for Pind in 1:P-1

            (iind,jind) = SP[Pind,:]

            temp += len_ij_U_ij_U_X_2U_X_U_svd[iind,jind,x4,y4]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_2U_X_U_svd[i1,j1,x4,y4,1:U1-1])

    end

    (a1,b1) = SU[i1,j1,U1,3:4]

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_U_X_2U_X_U_svd[a1,b1,i2ind,j2ind,x4,y4]

        end

    end

    if Z2 > 1

        temp += sum(len_U_X_2U_X_U_svd[a1,b1,i2,j2,x4,y4,1:Z2-1])

    end

    (x2,y2) = SU[i2,j2,Z2,3:4]

    if X2 > 1

        temp += sum(len_X_2U_X_U_svd[a1,b1,x2,y2,x4,y4,1:X2-1])

    end

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    if U4 > 1

        temp += sum(len_2U_X_U_svd[a2,b2,x4,y4,1:U4-1])

    end

    (a3,b3) = SU[a2,b2,U4,3:4]

    if U5 > 1

        temp += sum(len_U_X_U_svd[a3,b3,x4,y4,1:U5-1])

    end

    (a4,b4) = SU[a3,b3,U5,3:4]

    if X4 > 1

        temp += sum(len_X_U[a4,b4,x4,y4,1:X4-1])

    end

    (a5,b5) = SX[a4,b4,x4,y4,X4,:]

    temp += U8

    return temp

end

function testing(x4::Int,y4::Int)

    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind2_svd[x4,y4])
    println(@time translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind2_vec_to_ind_naive(test1,x4,y4))
    (i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81) = translate_10p_SVD_ind2_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind2_ind_to_vec(i11,j11,U11,i21,j21,Z21,X21,U41,U51,X41,U81,x4,y4) - test1)

end

#testing(1,1)
#testing(1,2)
#testing(2,1)
#testing(2,2)


function translate_10p_SVD_ind1_vec_to_ind(index::Int,x4::Int,y4::Int)

    #the index is simply an U

    Z4 = index

    return Z4

end

function translate_10p_SVD_ind1_vec_to_ind_naive(index::Int,x4::Int,y4::Int)

    Z40 = 0

    count = 0

    for Z4 in 1:RangeU[x4,y4]

        (a1,b1) = SU[x4,y4,Z4,3:4]

        count += 1

        if count == index

            Z40 = Z4


        end

    end

    return Z40

end


function translate_10p_SVD_ind1_ind_to_vec(Z4::Int,x4::Int,y4::Int)

    return Z4

end

function testing_2(x4::Int,y4::Int)

    test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    Z41 = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    Z41 = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    Z41 = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    Z41 = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,x4,y4) - test1)


    test1 = rand(1:configs_10p_ind1_svd[x4,y4])
    println(@time translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4))
    println(@time translate_10p_SVD_ind1_vec_to_ind_naive(test1,x4,y4))
    Z41 = translate_10p_SVD_ind1_vec_to_ind(test1,x4,y4)
    println(@time translate_10p_SVD_ind1_ind_to_vec(Z41,x4,y4) - test1)

end

#testing_2(1,1)
#testing_2(1,2)
#testing_2(2,1)
#testing_2(2,2)
