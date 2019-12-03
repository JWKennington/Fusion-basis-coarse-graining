function translate_8p_SVD_ind2_vec_to_ind(index::Int,x3::Int,y3::Int)

    #u2,K2,U1,s1,K1,U2,E1p,X1,X2,U4,E

    P2 = 0
    U1 = 0
    P1 = 0
    U2 = 0
    X1 = 0
    X2 = 0
    U4 = 0

    temp = index

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_U_ij_U_X_X_U_svd[i2,j2,x3,y3]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_U_ij_U_X_X_U_svd[i2,j2,x3,y3]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_X_U_svd[i2,j2,x3,y3,U1]

    end

    temp += len_U_ij_U_X_X_U_svd[i2,j2,x3,y3,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        P1 += 1

        (i1,j1) = SP[P1,:]

        temp -= len_ij_U_X_X_U_svd[a1,b1,i1,j1,x3,y3]

    end

    (i1,j1) = SP[P1,:]

    temp += len_ij_U_X_X_U_svd[a1,b1,i1,j1,x3,y3]

    while temp > 0

        U2 += 1

        temp -= len_U_X_X_U_svd[a1,b1,i1,j1,x3,y3,U2]

    end

    temp += len_U_X_X_U_svd[a1,b1,i1,j1,x3,y3,U2]

    (x1,y1) = SU[i1,j1,U2,3:4]

    while temp > 0

        X1 += 1

        temp -= len_X_X_U_svd[a1,b1,x1,y1,x3,y3,X1]

    end

    temp += len_X_X_U_svd[a1,b1,x1,y1,x3,y3,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        X2 += 1

        temp -= len_X_U[a2,b2,x3,y3,X2]

    end

    temp += len_X_U[a2,b2,x3,y3,X2]

    (a3,b3) = SX[a2,b2,x3,y3,X2,:]

    U4 = temp

    return (i2,j2,U1,i1,j1,U2,X1,X2,U4)

end



function translate_8p_SVD_ind2_vec_to_ind_naive(index::Int,x3::Int,y3::Int)

    i20 = 0
    j20 = 0
    U10 = 0
    i10 = 0
    j10 = 0
    U20 = 0
    X10 = 0
    X20 = 0
    U40 = 0

    count = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for i1 in 1:k+1, j1 in 1:k+1

                for U2 in 1:RangeU[i1,j1]

                    (k1,l1,x1,y1) = SU[i1,j1,U2,:]

                    for X1 in 1:RangeX[M1,N1,x1,y1]

                        (M5,N5) = SX[M1,N1,x1,y1,X1,:]

                        for X2 in 1:RangeX[M5,N5,x3,y3]

                            (M3,N3) = SX[M5,N5,x3,y3,X2,:]

                            for U4 in 1:RangeU[M3,N3]

                                (k6,l6,x4,y4) = SU[M3,N3,U4,:]

                                count += 1

                                if count == index

                                    (i20,j20,U10,i10,j10,U20,X10,X20,U40) =
                                        (i2,j2,U1,i1,j1,U2,X1,X2,U4)

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,j20,U10,i10,j10,U20,X10,X20,U40)

end



function translate_8p_SVD_ind2_ind_to_vec(i2::Int,j2::Int,U1::Int,
    i1::Int,j1::Int,U2::Int,X1::Int,
    X2::Int,U4::Int,
    x3::Int,y3::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_U_ij_U_X_X_U_svd[i2ind,j2ind,x3,y3]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_X_U_svd[i2,j2,x3,y3,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    P1 = Pf[i1,j1]

    if P1 > 1

        for P1ind in 1:P1-1

            (i1ind,j1ind) = SP[P1ind,:]

            temp += len_ij_U_X_X_U_svd[a1,b1,i1ind,j1ind,x3,y3]

        end

    end

    if U2 > 1

        temp += sum(len_U_X_X_U_svd[a1,b1,i1,j1,x3,y3,1:U2-1])

    end

    (x1,y1) = SU[i1,j1,U2,3:4]

    if X1 > 1

        temp += sum(len_X_X_U_svd[a1,b1,x1,y1,x3,y3,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if X2 > 1

        temp += sum(len_X_U[a2,b2,x3,y3,1:X2-1])

    end

    (a3,b3) = SX[a2,b2,x3,y3,X2,:]

    temp += U4

    return temp

end

function testing_8p(x3::Int,y3::Int)

    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i21,j21,U11,i11,j11,U21,X11,X21,U41) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i21,j21,U11,i11,j11,U21,X11,X21,U41,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i21,j21,U11,i11,j11,U21,X11,X21,U41) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i21,j21,U11,i11,j11,U21,X11,X21,U41,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i21,j21,U11,i11,j11,U21,X11,X21,U41) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i21,j21,U11,i11,j11,U21,X11,X21,U41,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i21,j21,U11,i11,j11,U21,X11,X21,U41) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i21,j21,U11,i11,j11,U21,X11,X21,U41,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind2_svd[x3,y3])
    println(@time translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind2_vec_to_ind_naive(test1,x3,y3))
    (i21,j21,U11,i11,j11,U21,X11,X21,U41) = translate_8p_SVD_ind2_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind2_ind_to_vec(i21,j21,U11,i11,j11,U21,X11,X21,U41,x3,y3) - test1)

end

#testing_8p(1,1)
#testing_8p(1,2)
#testing_8p(2,1)
#testing_8p(2,2)


function translate_8p_SVD_ind1_vec_to_ind(index::Int,x3::Int,y3::Int)

    #Z3 = 0
    #E3 = 0

    Z3 = index

    return Z3

end

function translate_8p_SVD_ind1_vec_to_ind_naive(index::Int,x3::Int,y3::Int)

    Z30 = 0

    count = 0

    for Z3 in 1:RangeU[x3,y3]

        (a1,b1) = SU[x3,y3,Z3,3:4]

        count += 1

        if count == index

            Z30 = Z3

        end

    end

    return Z30

end


function translate_8p_SVD_ind1_ind_to_vec(Z3::Int,x3::Int,y3::Int)

    temp = Z3

    return temp

end

function testing_8p_2(x3::Int,y3::Int)

    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
    Z31 = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
    Z31 = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
    Z31 = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
    Z31 = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,x3,y3) - test1)


    test1 = rand(1:configs_8p_ind1_svd[x3,y3])
    println(@time translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3))
    println(@time translate_8p_SVD_ind1_vec_to_ind_naive(test1,x3,y3))
    Z31 = translate_8p_SVD_ind1_vec_to_ind(test1,x3,y3)
    println(@time translate_8p_SVD_ind1_ind_to_vec(Z31,x3,y3) - test1)

end

#testing_2(1,1)
#testing_2(1,2)
#testing_2(2,1)
#testing_2(2,2)
