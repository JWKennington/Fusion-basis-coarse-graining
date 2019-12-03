#Define functions to translate between 8p states, indices and vectorized

function translate_8p_vec_to_ind_1(index::Int)
    #s1,K1,U1,U2,U3,U4,U5,U6,E

    temp = index

    P1 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    U5 = 0
    U6 = 0

    while temp > 0

        P1 += 1

        (i1,j1) = SP[P1,:]

        temp -= len_ij_6U[i1,j1]

    end

    (i1,j1) = SP[P1,:]

    temp += len_ij_6U[i1,j1]

    while temp > 0

        U1 += 1

        temp -= len_6U[i1,j1,U1]

    end

    temp += len_6U[i1,j1,U1]

    (a1,b1) = SU[i1,j1,U1,3:4]

    while temp > 0

        U2 += 1

        temp -= len_5U[a1,b1,U2]

    end

    temp += len_5U[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,3:4]

    while temp > 0

        U3 += 1

        temp -= len_4U[a2,b2,U3]

    end

    temp += len_4U[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,3:4]

    while temp > 0

        U4 += 1

        temp -= len_3U[a3,b3,U4]

    end

    temp += len_3U[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,3:4]

    while temp > 0

        U5 += 1

        temp -= len_2U[a4,b4,U5]

    end

    temp += len_2U[a4,b4,U5]

    (a5,b5) = SU[a4,b4,U5,3:4]

    U6 = temp

    return (i1,j1,U1,U2,U3,U4,U5,U6)

end


function translate_8p_vec_to_ind_1_naive(index::Int)
    #s1,K1,U1,U2,U3,U4,U5,U6,E

    count = 0

    i10 = 0
    j10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    U50 = 0
    U60 = 0

    for i1 in 1:k+1, j1 in 1:k+1

        for U1 in 1:RangeU[i1,j1]

            (i3,j3,a1,b1) = SU[i1,j1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (x2,y2,a3,b3) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a3,b3]

                    (i5,j5,M2,N2) = SU[a3,b3,U3,:]

                    for U4 in 1:RangeU[M2,N2]

                        (k6,l6,c2,d2) = SU[M2,N2,U4,:]

                        for U5 in 1:RangeU[c2,d2]

                            (x4,y4,c1,d1) = SU[c2,d2,U5,:]

                            for U6 in 1:RangeU[c1,d1]

                                (k3,l3,k1,l1) = SU[c1,d1,U6,:]

                                count += 1

                                if count == index

                                    (i10,j10,U10,U20,U30,U40,U50,U60) =
                                        (i1,j1,U1,U2,U3,U4,U5,U6)

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i10,j10,U10,U20,U30,U40,U50,U60)

end


function translate_8p_ind_to_vec_1(i1::Int,j1::Int,U1::Int,U2::Int,
    U3::Int,U4::Int,U5::Int,U6::Int)

    temp = 0

    P1 = Pf[i1,j1]

    if P1 > 1

        for P1ind in 1:P1-1

            (i1ind,j1ind) = SP[P1ind,:]

            temp += len_ij_6U[i1ind,j1ind]

        end

    end

    if U1 > 1

        temp += sum(len_6U[i1,j1,1:U1-1])

    end

    (a1,b1) = SU[i1,j1,U1,3:4]

    if U2 >1

        temp += sum(len_5U[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,3:4]

    if U3 > 1

        temp += sum(len_4U[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,3:4]

    if U4 > 1

        temp += sum(len_3U[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,3:4]

    if U5 > 1

        temp += sum(len_2U[a4,b4,1:U5-1])

    end

    (a5,b5) = SU[a4,b4,U5,3:4]

    temp += U6

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_1(test1))
#println(translate_8p_vec_to_ind_1_naive(test1))
#(i11,j11,U11,U21,U31,U41,U51,U61) = translate_8p_vec_to_ind_1(test1)
#println(translate_8p_ind_to_vec_1(i11,j11,U11,U21,U31,U41,U51,U61) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_1(test1))
#println(translate_8p_vec_to_ind_1_naive(test1))
#(i11,j11,U11,U21,U31,U41,U51,U61) = translate_8p_vec_to_ind_1(test1)
#println(translate_8p_ind_to_vec_1(i11,j11,U11,U21,U31,U41,U51,U61) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_1(test1))
#println(translate_8p_vec_to_ind_1_naive(test1))
#(i11,j11,U11,U21,U31,U41,U51,U61) = translate_8p_vec_to_ind_1(test1)
#println(translate_8p_ind_to_vec_1(i11,j11,U11,U21,U31,U41,U51,U61) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_1(test1))
#println(translate_8p_vec_to_ind_1_naive(test1))
#(i11,j11,U11,U21,U31,U41,U51,U61) = translate_8p_vec_to_ind_1(test1)
#println(translate_8p_ind_to_vec_1(i11,j11,U11,U21,U31,U41,U51,U61) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_1(test1))
#println(translate_8p_vec_to_ind_1_naive(test1))
#(i11,j11,U11,U21,U31,U41,U51,U61) = translate_8p_vec_to_ind_1(test1)
#println(translate_8p_ind_to_vec_1(i11,j11,U11,U21,U31,U41,U51,U61) - test1)


#Transformation 1 tested

function translate_8p_vec_to_ind_2(index::Int)
    #u2,K2,U1,s1,K1,U2,X1,U3,U4,U5,E

    temp = index

    P2 = 0
    U1 = 0
    P1 = 0
    U2 = 0
    X1 = 0
    U3 = 0
    U4 = 0
    U5 = 0

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_U_ij_U_X_3U[i2,j2]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_U_ij_U_X_3U[i2,j2]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_3U[i2,j2,U1]

    end

    temp += len_U_ij_U_X_3U[i2,j2,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        P1 += 1

        (i1,j1) = SP[P1,:]

        temp -= len_ij_U_X_3U[a1,b1,i1,j1]

    end

    (i1,j1) = SP[P1,:]

    temp += len_ij_U_X_3U[a1,b1,i1,j1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_3U[a1,b1,i1,j1,U2]

    end

    temp += len_U_X_3U[a1,b1,i1,j1,U2]

    (x1,y1) = SU[i1,j1,U2,3:4]

    while temp > 0

        X1 += 1

        temp -= len_X_3U[a1,b1,x1,y1,X1]

    end

    temp += len_X_3U[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        U3 += 1

        temp -= len_3U[a2,b2,U3]

    end

    temp += len_3U[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,3:4]

    while temp > 0

        U4 += 1

        temp -= len_2U[a3,b3,U4]

    end

    temp += len_2U[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,3:4]

    U5 = temp

    return (i2,j2,U1,i1,j1,U2,X1,U3,U4,U5)

end


function translate_8p_vec_to_ind_2_naive(index::Int)
    #u2,K2,U1,s1,K1,U2,X1,U3,U4,U5,E

    count = 0

    i20 = 0
    j20 = 0
    U10 = 0
    i10 = 0
    j10 = 0
    U20 = 0
    X10 = 0
    U30 = 0
    U40 = 0
    U50 = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for i1 in 1:k+1, j1 in 1:k+1

                for U2 in 1:RangeU[i1,j1]

                    (i3,j3,a1,b1) = SU[i1,j1,U2,:]

                    for X1 in 1:RangeX[M1,N1,a1,b1]

                        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

                        for U3 in 1:RangeU[M2,N2]

                            (k6,l6,c2,d2) = SU[M2,N2,U3,:]

                            for U4 in 1:RangeU[c2,d2]

                                (i4,j4,c1,d1) = SU[c2,d2,U4,:]

                                for U5 in 1:RangeU[c1,d1]

                                    (k3,l3,k1,l1) = SU[c1,d1,U5,:]

                                    count += 1

                                    if count == index

                                        (i20,j20,U10,i10,j10,U20,X10,U30,U40,U50) =
                                            (i2,j2,U1,i1,j1,U2,X1,U3,U4,U5)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,j20,U10,i10,j10,U20,X10,U30,U40,U50)

end


function translate_8p_ind_to_vec_2(i2::Int,j2::Int,U1::Int,
    i1::Int,j1::Int,U2::Int,X1::Int,
    U3::Int,U4::Int,U5::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_U_ij_U_X_3U[i2ind,j2ind]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_3U[i2,j2,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    P1 = Pf[i1,j1]

    if P1 > 1

        for P1ind in 1:P1-1

            (i1ind,j1ind) = SP[P1ind,:]

            temp += len_ij_U_X_3U[a1,b1,i1ind,j1ind]

        end

    end

    if U2 > 1

        temp += sum(len_U_X_3U[a1,b1,i1,j1,1:U2-1])

    end

    (x1,y1) = SU[i1,j1,U2,3:4]

    if X1 > 1

        temp += sum(len_X_3U[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if U3 > 1

        temp += sum(len_3U[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,3:4]

    if U4 > 1

        temp += sum(len_2U[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,3:4]

    temp += U5

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_2(test1))
#println(translate_8p_vec_to_ind_2_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) = translate_8p_vec_to_ind_2(test1)
#println(translate_8p_ind_to_vec_2(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_2(test1))
#println(translate_8p_vec_to_ind_2_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) = translate_8p_vec_to_ind_2(test1)
#println(translate_8p_ind_to_vec_2(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_2(test1))
#println(translate_8p_vec_to_ind_2_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) = translate_8p_vec_to_ind_2(test1)
#println(translate_8p_ind_to_vec_2(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_2(test1))
#println(translate_8p_vec_to_ind_2_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) = translate_8p_vec_to_ind_2(test1)
#println(translate_8p_ind_to_vec_2(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_2(test1))
#println(translate_8p_vec_to_ind_2_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) = translate_8p_vec_to_ind_2(test1)
#println(translate_8p_ind_to_vec_2(i21,j21,U11,i11,j11,U21,X11,U31,U41,U51) - test1)


#Translations for trafo_2


function translate_8p_vec_to_ind_3(index::Int)
    #u2,K2,U1,s1,K1,U2,X1,t1,L1,U3,X2,U4,E

    temp = index

    P2 = 0
    U1 = 0
    P1 = 0
    U2 = 0
    X1 = 0
    Q1 = 0
    U3 = 0
    X2 = 0
    U4 = 0

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_U_ij_U_X_ij_U_X_U[i2,j2]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_U_ij_U_X_ij_U_X_U[i2,j2]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_ij_U_X_U[i2,j2,U1]

    end

    temp += len_U_ij_U_X_ij_U_X_U[i2,j2,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        P1 += 1

        (i1,j1) = SP[P1,:]

        temp -= len_ij_U_X_ij_U_X_U[a1,b1,i1,j1]

    end

    (i1,j1) = SP[P1,:]

    temp += len_ij_U_X_ij_U_X_U[a1,b1,i1,j1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_ij_U_X_U[a1,b1,i1,j1,U2]

    end

    temp += len_U_X_ij_U_X_U[a1,b1,i1,j1,U2]

    (x1,y1) = SU[i1,j1,U2,3:4]

    while temp > 0

        X1 += 1

        temp -= len_X_ij_U_X_U[a1,b1,x1,y1,X1]

    end

    temp += len_X_ij_U_X_U[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        Q1 += 1

        (k1,l1) = SP[Q1,:]

        temp -= len_ij_U_X_U[a2,b2,k1,l1]

    end

    (k1,l1) = SP[Q1,:]

    temp += len_ij_U_X_U[a2,b2,k1,l1]

    while temp > 0

        U3 += 1

        temp -= len_U_X_U[a2,b2,k1,l1,U3]

    end

    temp += len_U_X_U[a2,b2,k1,l1,U3]

    (x2,y2) = SU[k1,l1,U3,3:4]

    while temp > 0

        X2 += 1

        temp -= len_X_U[a2,b2,x2,y2,X2]

    end

    temp += len_X_U[a2,b2,x2,y2,X2]

    (a3,b3) = SX[a2,b2,x2,y2,X2,:]

    U4 = temp

    return (i2,j2,U1,i1,j1,U2,X1,k1,l1,U3,X2,U4)

end


function translate_8p_vec_to_ind_3_naive(index::Int)
    #u2,K2,U1,s1,K1,U2,X1,t1,L1,U3,X2,U4,E

    count = 0

    i20 = 0
    j20 = 0
    U10 = 0
    i10 = 0
    j10 = 0
    U20 = 0
    X10 = 0
    k10 = 0
    l10 = 0
    U30 = 0
    X20 = 0
    U40 = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for i1 in 1:k+1, j1 in 1:k+1

                for U2 in 1:RangeU[i1,j1]

                    (i3,j3,a1,b1) = SU[i1,j1,U2,:]

                    for X1 in 1:RangeX[M1,N1,a1,b1]

                        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

                        for k1 in 1:k+1, l1 in 1:k+1

                            for U3 in 1:RangeU[k1,l1]

                                (k3,l3,c1,d1) = SU[k1,l1,U3,:]

                                for X2 in 1:RangeX[M2,N2,c1,d1]

                                    (M3,N3) = SX[M2,N2,c1,d1,X2,:]

                                    for U4 in 1:RangeU[M3,N3]

                                        (k6,l6,x4,y4) = SU[M3,N3,U4,:]

                                        count += 1

                                        if count == index

                                            (i20,j20,U10,i10,j10,U20,X10,k10,l10,U30,X20,U40) =
                                                (i2,j2,U1,i1,j1,U2,X1,k1,l1,U3,X2,U4)

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

    return (i20,j20,U10,i10,j10,U20,X10,k10,l10,U30,X20,U40)

end


function translate_8p_ind_to_vec_3(i2::Int,j2::Int,U1::Int,
    i1::Int,j1::Int,U2::Int,X1::Int,
    k1::Int,l1::Int,U3::Int,X2::Int,
    U4::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_U_ij_U_X_ij_U_X_U[i2ind,j2ind]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_ij_U_X_U[i2,j2,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    P1 = Pf[i1,j1]

    if P1 > 1

        for P1ind in 1:P1-1

            (i1ind,j1ind) = SP[P1ind,:]

            temp += len_ij_U_X_ij_U_X_U[a1,b1,i1ind,j1ind]

        end

    end

    if U2 > 1

        temp += sum(len_U_X_ij_U_X_U[a1,b1,i1,j1,1:U2-1])

    end

    (x1,y1) = SU[i1,j1,U2,3:4]

    if X1 > 1

        temp += sum(len_X_ij_U_X_U[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    Q1 = Pf[k1,l1]

    if Q1 > 1

        for Q1ind in 1:Q1-1

            (k1ind,l1ind) = SP[Q1ind,:]

            temp += len_ij_U_X_U[a2,b2,k1ind,l1ind]

        end

    end

    if U3 > 1

        temp += sum(len_U_X_U[a2,b2,k1,l1,1:U3-1])

    end

    (x2,y2) = SU[k1,l1,U3,3:4]

    if X2 > 1

        temp += sum(len_X_U[a2,b2,x2,y2,1:X2-1])

    end

    (a3,b3) = SX[a2,b2,x2,y2,X2,:]

    temp += U4

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_3(test1))
#println(translate_8p_vec_to_ind_3_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) = translate_8p_vec_to_ind_3(test1)
#println(translate_8p_ind_to_vec_3(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_3(test1))
#println(translate_8p_vec_to_ind_3_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) = translate_8p_vec_to_ind_3(test1)
#println(translate_8p_ind_to_vec_3(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_3(test1))
#println(translate_8p_vec_to_ind_3_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) = translate_8p_vec_to_ind_3(test1)
#println(translate_8p_ind_to_vec_3(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_3(test1))
#println(translate_8p_vec_to_ind_3_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) = translate_8p_vec_to_ind_3(test1)
#println(translate_8p_ind_to_vec_3(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_3(test1))
#println(translate_8p_vec_to_ind_3_naive(test1))
#(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) = translate_8p_vec_to_ind_3(test1)
#println(translate_8p_ind_to_vec_3(i21,j21,U11,i11,j11,U21,X11,k11,l11,U31,X21,U41) - test1)


#Functions for translating for trafo_3 and trafo_6


function translate_8p_vec_to_ind_4(index::Int)
    #u2,K2,U1,U2,U3,t1,L1,U4,X2,U5,E

    temp = index

    P2 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    Q1 = 0
    U4 = 0
    X2 = 0
    U5 = 0

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_3U_ij_U_X_U[i2,j2]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_3U_ij_U_X_U[i2,j2]

    while temp > 0

        U1 += 1

        temp -= len_3U_ij_U_X_U[i2,j2,U1]

    end

    temp += len_3U_ij_U_X_U[i2,j2,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        U2 += 1

        temp -= len_2U_ij_U_X_U[a1,b1,U2]

    end

    temp += len_2U_ij_U_X_U[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,3:4]

    while temp > 0

        U3 += 1

        temp -= len_U_ij_U_X_U[a2,b2,U3]

    end

    temp += len_U_ij_U_X_U[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,3:4]

    while temp > 0

        Q1 += 1

        (k1,l1) = SP[Q1,:]

        temp -= len_ij_U_X_U[a3,b3,k1,l1]

    end

    (k1,l1) = SP[Q1,:]

    temp += len_ij_U_X_U[a3,b3,k1,l1]

    while temp > 0

        U4 += 1

        temp -= len_U_X_U[a3,b3,k1,l1,U4]

    end

    temp += len_U_X_U[a3,b3,k1,l1,U4]

    (x2,y2) = SU[k1,l1,U4,3:4]

    while temp > 0

        X2 += 1

        temp -= len_X_U[a3,b3,x2,y2,X2]

    end

    temp += len_X_U[a3,b3,x2,y2,X2]

    (a4,b4) = SX[a3,b3,x2,y2,X2,:]

    U5 = temp

    return (i2,j2,U1,U2,U3,k1,l1,U4,X2,U5)

end


function translate_8p_vec_to_ind_4_naive(index::Int)
    #u2,K2,U1,U2,U3,t1,L1,U4,X2,U5,E

    count = 0

    i20 = 0
    j20 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    k10 = 0
    l10 = 0
    U40 = 0
    X20 = 0
    U50 = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for U2 in 1:RangeU[M1,N1]

                (i1,j1,M4,N4) = SU[M1,N1,U2,:]

                for U3 in 1:RangeU[M4,N4]

                    (i3,j3,M2,N2) = SU[M4,N4,U3,:]

                    for k1 in 1:k+1, l1 in 1:k+1

                        for U4 in 1:RangeU[k1,l1]

                            (k3,l3,c1,d1) = SU[k1,l1,U4,:]

                            for X2 in 1:RangeX[M2,N2,c1,d1]

                                (M3,N3) = SX[M2,N2,c1,d1,X2,:]

                                for U5 in 1:RangeU[M3,N3]

                                    (k6,l6,x4,y4) = SU[M3,N3,U5,:]

                                    count += 1

                                    if count == index

                                        (i20,j20,U10,U20,U30,k10,l10,U40,X20,U50) =
                                            (i2,j2,U1,U2,U3,k1,l1,U4,X2,U5)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,j20,U10,U20,U30,k10,l10,U40,X20,U50)

end


function translate_8p_ind_to_vec_4(i2::Int,j2::Int,U1::Int,U2::Int,U3::Int,
    k1::Int,l1::Int,U4::Int,X2::Int,
    U5::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_3U_ij_U_X_U[i2ind,j2ind]

        end

    end

    if U1 > 1

        temp += sum(len_3U_ij_U_X_U[i2,j2,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    if U2 > 1

        temp += sum(len_2U_ij_U_X_U[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,3:4]

    if U3 > 1

        temp += sum(len_U_ij_U_X_U[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,3:4]

    Q1 = Pf[k1,l1]

    if Q1 > 1

        for Q1ind in 1:Q1-1

            (k1ind,l1ind) = SP[Q1ind,:]

            temp += len_ij_U_X_U[a3,b3,k1ind,l1ind]

        end

    end

    if U4 > 1

        temp += sum(len_U_X_U[a3,b3,k1,l1,1:U4-1])

    end

    (x1,y1) = SU[k1,l1,U4,3:4]

    if X2 > 1

        temp += sum(len_X_U[a3,b3,x1,y1,1:X2-1])

    end

    (a4,b4) = SX[a3,b3,x1,y1,X2,:]

    temp += U5

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_4(test1))
#println(translate_8p_vec_to_ind_4_naive(test1))
#(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) = translate_8p_vec_to_ind_4(test1)
#println(translate_8p_ind_to_vec_4(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_4(test1))
#println(translate_8p_vec_to_ind_4_naive(test1))
#(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) = translate_8p_vec_to_ind_4(test1)
#println(translate_8p_ind_to_vec_4(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_4(test1))
#println(translate_8p_vec_to_ind_4_naive(test1))
#(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) = translate_8p_vec_to_ind_4(test1)
#println(translate_8p_ind_to_vec_4(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_4(test1))
#println(translate_8p_vec_to_ind_4_naive(test1))
#(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) = translate_8p_vec_to_ind_4(test1)
#println(translate_8p_ind_to_vec_4(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_4(test1))
#println(translate_8p_vec_to_ind_4_naive(test1))
#(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) = translate_8p_vec_to_ind_4(test1)
#println(translate_8p_ind_to_vec_4(i21,j21,U11,U21,U31,k11,l11,U41,X21,U51) - test1)


# translate functions for trafo_4 and trafo_5

function translate_8p_vec_to_ind_5(index::Int)
    #u2,K2,U1,U2,t1,L1,U3,U4,X2,U5,E

    temp = index

    P2 = 0
    U1 = 0
    U2 = 0
    Q1 = 0
    U3 = 0
    U4 = 0
    X2 = 0
    U5 = 0

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_2U_ij_2U_X_U[i2,j2]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_2U_ij_2U_X_U[i2,j2]

    while temp > 0

        U1 += 1

        temp -= len_2U_ij_2U_X_U[i2,j2,U1]

    end

    temp += len_2U_ij_2U_X_U[i2,j2,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        U2 += 1

        temp -= len_U_ij_2U_X_U[a1,b1,U2]

    end

    temp += len_U_ij_2U_X_U[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,3:4]

    while temp > 0

        Q1 += 1

        (k1,l1) = SP[Q1,:]

        temp -= len_ij_2U_X_U[a2,b2,k1,l1]

    end

    (k1,l1) = SP[Q1,:]

    temp += len_ij_2U_X_U[a2,b2,k1,l1]

    while temp > 0

        U3 += 1

        temp -= len_2U_X_U[a2,b2,k1,l1,U3]

    end

    temp += len_2U_X_U[a2,b2,k1,l1,U3]

    (a3,b3) = SU[k1,l1,U3,3:4]

    while temp > 0

        U4 += 1

        temp -= len_U_X_U[a2,b2,a3,b3,U4]

    end

    temp += len_U_X_U[a2,b2,a3,b3,U4]

    (x2,y2) = SU[a3,b3,U4,3:4]

    while temp > 0

        X2 += 1

        temp -= len_X_U[a2,b2,x2,y2,X2]

    end

    temp += len_X_U[a2,b2,x2,y2,X2]

    (a4,b4) = SX[a2,b2,x2,y2,X2,:]

    U5 = temp

    return (i2,j2,U1,U2,k1,l1,U3,U4,X2,U5)

end


function translate_8p_vec_to_ind_5_naive(index::Int)
    #u2,K2,U1,U2,t1,L1,U3,U4,X2,U5,E

    count = 0

    i20 = 0
    j20 = 0
    U10 = 0
    U20 = 0
    k10 = 0
    l10 = 0
    U30 = 0
    U40 = 0
    X20 = 0
    U50 = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for U2 in 1:RangeU[M1,N1]

                (i1,j1,M4,N4) = SU[M1,N1,U2,:]

                for k1 in 1:k+1, l1 in 1:k+1

                    for U3 in 1:RangeU[k1,l1]

                        (k3,l3,c1,d1) = SU[k1,l1,U3,:]

                        for U4 in 1:RangeU[c1,d1]

                            (i3,j3,c2,d2) = SU[c1,d1,U4,:]

                            for X2 in 1:RangeX[M4,N4,c2,d2]

                                (M3,N3) = SX[M4,N4,c2,d2,X2,:]

                                for U5 in 1:RangeU[M3,N3]

                                    (k6,l6,x4,y4) = SU[M3,N3,U5,:]

                                    count += 1

                                    if count == index

                                        (i20,j20,U10,U20,k10,l10,U30,U40,X20,U50) =
                                            (i2,j2,U1,U2,k1,l1,U3,U4,X2,U5)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,j20,U10,U20,k10,l10,U30,U40,X20,U50)

end


function translate_8p_ind_to_vec_5(i2::Int,j2::Int,U1::Int,U2::Int,
    k1::Int,l1::Int,U3::Int,U4::Int,X2::Int,U5::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_2U_ij_2U_X_U[i2ind,j2ind]

        end

    end

    if U1 > 1

        temp += sum(len_2U_ij_2U_X_U[i2,j2,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    if U2 > 1

        temp += sum(len_U_ij_2U_X_U[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,3:4]

    Q1 = Pf[k1,l1]

    if Q1 > 1

        for Q1ind in 1:Q1-1

            (k1ind,l1ind) = SP[Q1ind,:]

            temp += len_ij_2U_X_U[a2,b2,k1ind,l1ind]

        end

    end

    if U3 > 1

        temp += sum(len_2U_X_U[a2,b2,k1,l1,1:U3-1])

    end

    (a3,b3) = SU[k1,l1,U3,3:4]

    if U4 > 1

        temp += sum(len_U_X_U[a2,b2,a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,3:4]

    if X2 > 1

        temp += sum(len_X_U[a2,b2,a4,b4,1:X2-1])

    end

    (a5,b5) = SX[a2,b2,a4,b4,X2,:]

    temp += U5

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i21,j21,U11,U21,k11,l11,U31,U41,X21,U51) - test1)


function translate_8p_vec_to_ind_bSVD(index::Int)

    #u2,K2,U1,s1,K1,U2,E1,X1,s3,K3,U3,E3,X2,U4,E

    temp = index

    P2 = 0
    U1 = 0
    P1 = 0
    U2 = 0
    X1 = 0
    P3 = 0
    U3 = 0
    X2 = 0
    U4 = 0

    while temp > 0

        P2 += 1

        (i2,j2) = SP[P2,:]

        temp -= len_ij_U_ij_U_X_ij_U_X_U[i2,j2]

    end

    (i2,j2) = SP[P2,:]

    temp += len_ij_U_ij_U_X_ij_U_X_U[i2,j2]

    while temp > 0

        U1 += 1

        temp -= len_U_ij_U_X_ij_U_X_U[i2,j2,U1]

    end

    temp += len_U_ij_U_X_ij_U_X_U[i2,j2,U1]

    (a1,b1) = SU[i2,j2,U1,3:4]

    while temp > 0

        P1 += 1

        (i1,j1) = SP[P1,:]

        temp -= len_ij_U_X_ij_U_X_U[a1,b1,i1,j1]

    end

    (i1,j1) = SP[P1,:]

    temp += len_ij_U_X_ij_U_X_U[a1,b1,i1,j1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_ij_U_X_U[a1,b1,i1,j1,U2]

    end

    temp += len_U_X_ij_U_X_U[a1,b1,i1,j1,U2]

    (x1,y1) = SU[i1,j1,U2,3:4]

    while temp > 0

        X1 += 1

        temp -= len_X_ij_U_X_U[a1,b1,x1,y1,X1]

    end

    temp += len_X_ij_U_X_U[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        P3 += 1

        (i3,j3) = SP[P3,:]

        temp -= len_ij_U_X_U[a2,b2,i3,j3]

    end

    (i3,j3) = SP[P3,:]

    temp += len_ij_U_X_U[a2,b2,i3,j3]

    while temp > 0

        U3 += 1

        temp -= len_U_X_U[a2,b2,i3,j3,U3]

    end

    temp += len_U_X_U[a2,b2,i3,j3,U3]

    (x3,y3) = SU[i3,j3,U3,3:4]

    while temp > 0

        X2 += 1

        temp -= len_X_U[a2,b2,x3,y3,X2]

    end

    temp += len_X_U[a2,b2,x3,y3,X2]

    (a3,b3) = SX[a2,b2,x3,y3,X2,:]

    U4 = temp

    return (i2,j2,U1,i1,j1,U2,X1,i3,j3,U3,X2,U4)

end


function translate_8p_vec_to_ind_naive_bSVD(index::Int)

    count = 0

    i20 = 0
    j20 = 0
    U10 = 0
    i10 = 0
    j10 = 0
    U20 = 0
    X10 = 0
    i30 = 0
    j30 = 0
    U30 = 0
    X20 = 0
    U40 = 0

    for i2 in 1:k+1, j2 in 1:k+1

        for U1 in 1:RangeU[i2,j2]

            (i5,j5,M1,N1) = SU[i2,j2,U1,:]

            for i1 in 1:k+1, j1 in 1:k+1

                for U2 in 1:RangeU[i1,j1]

                    (k1,l1,x1,y1) = SU[i1,j1,U2,:]

                    for X1 in 1:RangeX[M1,N1,x1,y1]

                        (M5,N5) = SX[M1,N1,x1,y1,X1,:]

                        for i3 in 1:k+1, j3 in 1:k+1

                            for U3 in 1:RangeU[i3,j3]

                                (k3,l3,x3,y3) = SU[i3,j3,U3,:]

                                for X2 in 1:RangeX[M5,N5,x3,y3]

                                    (M3,N3) = SX[M5,N5,x3,y3,X2,:]

                                    for U4 in 1:RangeU[M3,N3]

                                        (k6,l6,x4,y4) = SU[M3,N3,U4,:]

                                        count += 1

                                        if count == index

                                            (i20,j20,U10,i10,j10,U20,X10,i30,j30,U30,X20,U40) =
                                                (i2,j2,U1,i1,j1,U2,X1,i3,j3,U3,X2,U4)

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

    return (i20,j20,U10,i10,j10,U20,X10,i30,j30,U30,X20,U40)

end

#Define the translation from set of indices to the vectorized index.

function translate_8p_ind_to_vec_bSVD(i2::Int,j2::Int,U1::Int,
    i1::Int,j1::Int,U2::Int,X1::Int,
    i3::Int,j3::Int,U3::Int,X2::Int,
    U4::Int)

    temp = 0

    P2 = Pf[i2,j2]

    if P2 > 1

        for P2ind in 1:P2-1

            (i2ind,j2ind) = SP[P2ind,:]

            temp += len_ij_U_ij_U_X_ij_U_X_U[i2ind,j2ind]

        end

    end

    if U1 > 1

        temp += sum(len_U_ij_U_X_ij_U_X_U[i2,j2,1:U1-1])

    end

    (a1,b1) = SU[i2,j2,U1,3:4]

    P1 = Pf[i1,j1]

    if P1 > 1

        for P1ind in 1:P1-1

            (i1ind,j1ind) = SP[P1ind,:]

            temp += len_ij_U_X_ij_U_X_U[a1,b1,i1ind,j1ind]

        end

    end

    if U2 > 1

        temp += sum(len_U_X_ij_U_X_U[a1,b1,i1,j1,1:U2-1])

    end

    (x1,y1) = SU[i1,j1,U2,3:4]

    if X1 > 1

        temp += sum(len_X_ij_U_X_U[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    P3 = Pf[i3,j3]

    if P3 > 1

        for P3ind in 1:P3-1

            (i3ind,j3ind) = SP[P3ind,:]

            temp += len_ij_U_X_U[a2,b2,i3ind,j3ind]

        end

    end

    if U3 > 1

        temp += sum(len_U_X_U[a2,b2,i3,j3,1:U3-1])

    end

    (x3,y3) = SU[i3,j3,U3,3:4]

    if X2 > 1

        temp += sum(len_X_U[a2,b2,x3,y3,1:X2-1])

    end

    (a3,b3) = SX[a2,b2,x3,y3,X2,:]

    temp += U4

    return temp

end


#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(i21,j21,U11,i11,j11,U21,X11,i31,j31,U31,X21,U41) - test1)
