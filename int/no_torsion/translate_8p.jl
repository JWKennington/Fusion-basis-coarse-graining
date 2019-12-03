#Define functions to translate between 8p states, indices and vectorized

function translate_8p_vec_to_ind_1(index::Int)
    #i1,U1,U2,U3,U4,U5,F

    temp = index

    i1 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    U5 = 0
    F = 0

    while temp > 0

        i1 += 1

        temp -= len_i_5U_F[i1]

    end

    temp += len_i_5U_F[i1]

    while temp > 0

        U1 += 1

        temp -= len_5U_F[i1,i1,U1]

    end

    temp += len_5U_F[i1,i1,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_4U_F[a1,b1,U2]

    end

    temp += len_4U_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_3U_F[a2,b2,U3]

    end

    temp += len_3U_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_2U_F[a3,b3,U4]

    end

    temp += len_2U_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    while temp > 0

        U5 += 1

        temp -= len_U_F[a4,b4,U5]

    end

    temp += len_U_F[a4,b4,U5]

    F = temp

    return (i1,U1,U2,U3,U4,U5,F)

end


function translate_8p_vec_to_ind_1_naive(index::Int)
    #s1,K1,U1,U2,U3,U4,U5,U6,E

    count = 0

    i10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    U50 = 0
    F0 = 0

    for i1 in 1:k+1

        for U1 in 1:RangeU[i1,i1]

            (i3,a1,b1) = SU[i1,i1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (x2,a3,b3) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a3,b3]

                    (i5,M2,N2) = SU[a3,b3,U3,:]

                    for U4 in 1:RangeU[M2,N2]

                        (k6,c2,d2) = SU[M2,N2,U4,:]

                        for U5 in 1:RangeU[c2,d2]

                            (X4,c1,d1) = SU[c2,d2,U5,:]

                            for F in 1:RangeF[c1,d1]

                                (k3,k1) = SF[c1,d1,F,:]

                                count += 1

                                if count == index

                                    (i10,U10,U20,U30,U40,U50,F0) =
                                        (i1,U1,U2,U3,U4,U5,F)

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i10,U10,U20,U30,U40,U50,F0)

end


function translate_8p_ind_to_vec_1(i1::Int,U1::Int,U2::Int,
    U3::Int,U4::Int,U5::Int,F::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_5U_F[1:i1-1])

    end

    if U1 > 1

        temp += sum(len_5U_F[i1,i1,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if U2 >1

        temp += sum(len_4U_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_3U_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_2U_F[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if U5 > 1

        temp += sum(len_U_F[a4,b4,1:U5-1])

    end

    (a5,b5) = SU[a4,b4,U5,2:3]

    temp += F

    return temp

end



#Transformation 1 tested

function translate_8p_vec_to_ind_2(index::Int)
    #x2,U1,i1,U2,X1,U3,U4,F

    temp = index

    i2 = 0
    U1 = 0
    i1 = 0
    U2 = 0
    X1 = 0
    U3 = 0
    U4 = 0
    F = 0

    while temp > 0

        i2 += 1

        temp -= len_i_U_i_U_X_2U_F[i2]

    end

    temp += len_i_U_i_U_X_2U_F[i2]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_2U_F[i2,i2,U1]

    end

    temp += len_U_i_U_X_2U_F[i2,i2,U1]

    (a1,b1) = SU[i2,i2,U1,2:3]

    while temp > 0

        i1 += 1

        temp -= len_i_U_X_2U_F[a1,b1,i1]

    end

    temp += len_i_U_X_2U_F[a1,b1,i1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_2U_F[a1,b1,i1,i1,U2]

    end

    temp += len_U_X_2U_F[a1,b1,i1,i1,U2]

    (x1,y1) = SU[i1,i1,U2,2:3]

    while temp > 0

        X1 += 1

        temp -= len_X_2U_F[a1,b1,x1,y1,X1]

    end

    temp += len_X_2U_F[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        U3 += 1

        temp -= len_2U_F[a2,b2,U3]

    end

    temp += len_2U_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_U_F[a3,b3,U4]

    end

    temp += len_U_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    F = temp

    return (i2,U1,i1,U2,X1,U3,U4,F)

end


function translate_8p_vec_to_ind_2_naive(index::Int)


    count = 0

    i20 = 0
    U10 = 0
    i10 = 0
    U20 = 0
    X10 = 0
    U30 = 0
    U40 = 0
    F0 = 0

    for i2 in 1:k+1

        for U1 in 1:RangeU[i2,i2]

            (i5,M1,N1) = SU[i2,i2,U1,:]

            for i1 in 1:k+1

                for U2 in 1:RangeU[i1,i1]

                    (i3,a1,b1) = SU[i1,i1,U2,:]

                    for X1 in 1:RangeX[M1,N1,a1,b1]

                        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

                        for U3 in 1:RangeU[M2,N2]

                            (k6,c2,d2) = SU[M2,N2,U3,:]

                            for U4 in 1:RangeU[c2,d2]

                                (x4,c1,d1) = SU[c2,d2,U4,:]

                                for F in 1:RangeF[c1,d1]

                                    (k3,k1) = SF[c1,d1,F,:]

                                    count += 1

                                    if count == index

                                        (i20,U10,i10,U20,X10,U30,U40,F0) =
                                            (i2,U1,i1,U2,X1,U3,U4,F)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,U10,i10,U20,X10,U30,U40,F0)

end


function translate_8p_ind_to_vec_2(i2::Int,U1::Int,
    i1::Int,U2::Int,X1::Int,
    U3::Int,U4::Int,F::Int)

    temp = 0

    if i2 > 1

        temp += sum(len_i_U_i_U_X_2U_F[1:i2-1])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_2U_F[i2,i2,1:U1-1])

    end

    (a1,b1) = SU[i2,i2,U1,2:3]

    if i1 >1

        temp += sum(len_i_U_X_2U_F[a1,b1,1:i1-1])

    end

    if U2 > 1

        temp += sum(len_U_X_2U_F[a1,b1,i1,i1,1:U2-1])

    end

    (x1,y1) = SU[i1,i1,U2,2:3]

    if X1 > 1

        temp += sum(len_X_2U_F[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if U3 > 1

        temp += sum(len_2U_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_U_F[a3,b3,1:U4-1])

    end

    temp += F

    return temp

end


#Translations for trafo_2


function translate_8p_vec_to_ind_3(index::Int)
    #u2,K2,U1,s1,K1,U2,X1,t1,L1,U3,X2,U4,E

    temp = index

    i2 = 0
    U1 = 0
    i1 = 0
    U2 = 0
    X1 = 0
    k1 = 0
    U3 = 0
    X2 = 0
    F = 0

    while temp > 0

        i2 += 1

        temp -= len_i_U_i_U_X_i_U_X_F[i2]

    end

    temp += len_i_U_i_U_X_i_U_X_F[i2]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_i_U_X_F[i2,i2,U1]

    end

    temp += len_U_i_U_X_i_U_X_F[i2,i2,U1]

    (a1,b1) = SU[i2,i2,U1,2:3]

    while temp > 0

        i1 += 1

        temp -= len_i_U_X_i_U_X_F[a1,b1,i1]

    end

    temp += len_i_U_X_i_U_X_F[a1,b1,i1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_i_U_X_F[a1,b1,i1,i1,U2]

    end

    temp += len_U_X_i_U_X_F[a1,b1,i1,i1,U2]

    (x1,y1) = SU[i1,i1,U2,2:3]

    while temp > 0

        X1 += 1

        temp -= len_X_i_U_X_F[a1,b1,x1,y1,X1]

    end

    temp += len_X_i_U_X_F[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        k1 += 1

        temp -= len_i_U_X_F[a2,b2,k1]

    end

    temp += len_i_U_X_F[a2,b2,k1]

    while temp > 0

        U3 += 1

        temp -= len_U_X_F[a2,b2,k1,k1,U3]

    end

    temp += len_U_X_F[a2,b2,k1,k1,U3]

    (x2,y2) = SU[k1,k1,U3,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_F[a2,b2,x2,y2,X2]

    end

    temp += len_X_F[a2,b2,x2,y2,X2]

    F = temp

    return (i2,U1,i1,U2,X1,k1,U3,X2,F)

end


function translate_8p_vec_to_ind_3_naive(index::Int)

    count = 0

    i20 = 0
    U10 = 0
    i10 = 0
    U20 = 0
    X10 = 0
    k10 = 0
    U30 = 0
    X20 = 0
    F0 = 0

    for i2 in 1:k+1

        for U1 in 1:RangeU[i2,i2]

            (i5,M1,N1) = SU[i2,i2,U1,:]

            for i1 in 1:k+1

                for U2 in 1:RangeU[i1,i1]

                    (i3,a1,b1) = SU[i1,i1,U2,:]

                    for X1 in 1:RangeX[M1,N1,a1,b1]

                        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

                        for k1 in 1:k+1

                            for U3 in 1:RangeU[k1,k1]

                                (k3,c1,d1) = SU[k1,k1,U3,:]

                                for X2 in 1:RangeX[M2,N2,c1,d1]

                                    (M3,N3) = SX[M2,N2,c1,d1,X2,:]

                                    for F in 1:RangeF[M3,N3]

                                        (k6,x4) = SF[M3,N3,F,:]

                                        count += 1

                                        if count == index

                                            (i20,U10,i10,U20,X10,k10,U30,X20,F0) =
                                                (i2,U1,i1,U2,X1,k1,U3,X2,F)

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

    return (i20,U10,i10,U20,X10,k10,U30,X20,F0)

end


function translate_8p_ind_to_vec_3(i2::Int,U1::Int,
    i1::Int,U2::Int,X1::Int,
    k1::Int,U3::Int,X2::Int,F::Int)

    temp = 0

    if i2 > 1

        temp += sum(len_i_U_i_U_X_i_U_X_F[1:i2-1])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_i_U_X_F[i2,i2,1:U1-1])

    end

    (a1,b1) = SU[i2,i2,U1,2:3]

    if i1 >1

        temp += sum(len_i_U_X_i_U_X_F[a1,b1,1:i1-1])

    end

    if U2 > 1

        temp += sum(len_U_X_i_U_X_F[a1,b1,i1,i1,1:U2-1])

    end

    (x1,y1) = SU[i1,i1,U2,2:3]

    if X1 > 1

        temp += sum(len_X_i_U_X_F[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if k1 > 1

        temp += sum(len_i_U_X_F[a2,b2,1:k1-1])

    end

    if U3 > 1

        temp += sum(len_U_X_F[a2,b2,k1,k1,1:U3-1])

    end

    (x2,y2) = SU[k1,k1,U3,2:3]

    if X2 > 1

        temp += sum(len_X_F[a2,b2,x2,y2,1:X2-1])

    end

    temp += F

    return temp

end



#Functions for translating for trafo_3 and trafo_6


function translate_8p_vec_to_ind_4(index::Int)

    temp = index

    i2 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    k1 = 0
    U4 = 0
    X2 = 0
    F = 0

    while temp > 0

        i2 += 1

        temp -= len_i_3U_i_U_X_F[i2]

    end

    temp += len_i_3U_i_U_X_F[i2]

    while temp > 0

        U1 += 1

        temp -= len_3U_i_U_X_F[i2,i2,U1]

    end

    temp += len_3U_i_U_X_F[i2,i2,U1]

    (a1,b1) = SU[i2,i2,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_2U_i_U_X_F[a1,b1,U2]

    end

    temp += len_2U_i_U_X_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_U_i_U_X_F[a2,b2,U3]

    end

    temp += len_U_i_U_X_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        k1 += 1

        temp -= len_i_U_X_F[a3,b3,k1]

    end

    temp += len_i_U_X_F[a3,b3,k1]

    while temp > 0

        U4 += 1

        temp -= len_U_X_F[a3,b3,k1,k1,U4]

    end

    temp += len_U_X_F[a3,b3,k1,k1,U4]

    (x2,y2) = SU[k1,k1,U4,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_F[a3,b3,x2,y2,X2]

    end

    temp += len_X_F[a3,b3,x2,y2,X2]

    (a4,b4) = SX[a3,b3,x2,y2,X2,:]

    F = temp

    return (i2,U1,U2,U3,k1,U4,X2,F)

end


function translate_8p_vec_to_ind_4_naive(index::Int)

    count = 0

    i20 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    k10 = 0
    U40 = 0
    X20 = 0
    F0 = 0

    for i2 in 1:k+1

        for U1 in 1:RangeU[i2,i2]

            (i5,M1,N1) = SU[i2,i2,U1,:]

            for U2 in 1:RangeU[M1,N1]

                (i1,M4,N4) = SU[M1,N1,U2,:]

                for U3 in 1:RangeU[M4,N4]

                    (i3,M2,N2) = SU[M4,N4,U3,:]

                    for k1 in 1:k+1

                        for U4 in 1:RangeU[k1,k1]

                            (k3,c1,d1) = SU[k1,k1,U4,:]

                            for X2 in 1:RangeX[M2,N2,c1,d1]

                                (M3,N3) = SX[M2,N2,c1,d1,X2,:]

                                for F in 1:RangeF[M3,N3]

                                    (k6,x4) = SF[M3,N3,F,:]

                                    count += 1

                                    if count == index

                                        (i20,U10,U20,U30,k10,U40,X20,F0) =
                                            (i2,U1,U2,U3,k1,U4,X2,F)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,U10,U20,U30,k10,U40,X20,F0)

end


function translate_8p_ind_to_vec_4(i2::Int,U1::Int,U2::Int,
    U3::Int,k1::Int,U4::Int,X2,F)

    temp = 0

    if i2 > 1

        temp += sum(len_i_3U_i_U_X_F[1:i2-1])

    end

    if U1 > 1

        temp += sum(len_3U_i_U_X_F[i2,i2,1:U1-1])

    end

    (a1,b1) = SU[i2,i2,U1,2:3]

    if U2 > 1

        temp += sum(len_2U_i_U_X_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_U_i_U_X_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if k1 > 1

        temp += sum(len_i_U_X_F[a3,b3,1:k1-1])

    end

    if U4 > 1

        temp += sum(len_U_X_F[a3,b3,k1,k1,1:U4-1])

    end

    (x1,y1) = SU[k1,k1,U4,2:3]

    if X2 > 1

        temp += sum(len_X_F[a3,b3,x1,y1,1:X2-1])

    end

    temp += F

    return temp

end



# translate functions for trafo_4 and trafo_5

function translate_8p_vec_to_ind_5(index::Int)

    temp = index

    i2 = 0
    U1 = 0
    U2 = 0
    k1 = 0
    U3 = 0
    U4 = 0
    X2 = 0
    F = 0

    while temp > 0

        i2 += 1

        temp -= len_i_2U_i_2U_X_F[i2]

    end

    temp += len_i_2U_i_2U_X_F[i2]

    while temp > 0

        U1 += 1

        temp -= len_2U_i_2U_X_F[i2,i2,U1]

    end

    temp += len_2U_i_2U_X_F[i2,i2,U1]

    (a1,b1) = SU[i2,i2,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_U_i_2U_X_F[a1,b1,U2]

    end

    temp += len_U_i_2U_X_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        k1 += 1

        temp -= len_i_2U_X_F[a2,b2,k1]

    end

    temp += len_i_2U_X_F[a2,b2,k1]

    while temp > 0

        U3 += 1

        temp -= len_2U_X_F[a2,b2,k1,k1,U3]

    end

    temp += len_2U_X_F[a2,b2,k1,k1,U3]

    (a3,b3) = SU[k1,k1,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_U_X_F[a2,b2,a3,b3,U4]

    end

    temp += len_U_X_F[a2,b2,a3,b3,U4]

    (x2,y2) = SU[a3,b3,U4,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_F[a2,b2,x2,y2,X2]

    end

    temp += len_X_F[a2,b2,x2,y2,X2]

    F = temp

    return (i2,U1,U2,k1,U3,U4,X2,F)

end


function translate_8p_vec_to_ind_5_naive(index::Int)

    count = 0

    i20 = 0
    U10 = 0
    U20 = 0
    k10 = 0
    U30 = 0
    U40 = 0
    X20 = 0
    F0 = 0

    for i2 in 1:k+1

        for U1 in 1:RangeU[i2,i2]

            (i5,M1,N1) = SU[i2,i2,U1,:]

            for U2 in 1:RangeU[M1,N1]

                (i1,M4,N4) = SU[M1,N1,U2,:]

                for k1 in 1:k+1

                    for U3 in 1:RangeU[k1,k1]

                        (k3,c1,d1) = SU[k1,k1,U3,:]

                        for U4 in 1:RangeU[c1,d1]

                            (i3,c2,d2) = SU[c1,d1,U4,:]

                            for X2 in 1:RangeX[M4,N4,c2,d2]

                                (M3,N3) = SX[M4,N4,c2,d2,X2,:]

                                for F in 1:RangeF[M3,N3]

                                    (k6,x4) = SF[M3,N3,F,:]

                                    count += 1

                                    if count == index

                                        (i20,U10,U20,k10,U30,U40,X20,F0) =
                                            (i2,U1,U2,k1,U3,U4,X2,F)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return (i20,U10,U20,k10,U30,U40,X20,F0)

end


function translate_8p_ind_to_vec_5(i2::Int,U1::Int,
    U2::Int,k1::Int,U3::Int,U4::Int,X2::Int,F::Int)

    temp = 0

    if i2 > 1

        temp += sum(len_i_2U_i_2U_X_F[1:i2-1])

    end

    if U1 > 1

        temp += sum(len_2U_i_2U_X_F[i2,i2,1:U1-1])

    end

    (a1,b1) = SU[i2,i2,U1,2:3]

    if U2 > 1

        temp += sum(len_U_i_2U_X_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if k1 > 1

        temp += sum(len_i_2U_X_F[a2,b2,1:k1-1])

    end

    if U3 > 1

        temp += sum(len_2U_X_F[a2,b2,k1,k1,1:U3-1])

    end

    (a3,b3) = SU[k1,k1,U3,2:3]

    if U4 > 1

        temp += sum(len_U_X_F[a2,b2,a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if X2 > 1

        temp += sum(len_X_F[a2,b2,a4,b4,1:X2-1])

    end

    temp += F

    return temp

end


#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i2,U1,U2,k1,U3,U4,X2,F) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i2,U1,U2,k1,U3,U4,X2,F) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i2,U1,U2,k1,U3,U4,X2,F) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i2,U1,U2,k1,U3,U4,X2,F) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i2,U1,U2,k1,U3,U4,X2,F) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i2,U1,U2,k1,U3,U4,X2,F) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i2,U1,U2,k1,U3,U4,X2,F) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i2,U1,U2,k1,U3,U4,X2,F) - test1)

#test1 = rand(1:configs_8p)
#println(translate_8p_vec_to_ind_5(test1))
#println(translate_8p_vec_to_ind_5_naive(test1))
#(i2,U1,U2,k1,U3,U4,X2,F) = translate_8p_vec_to_ind_5(test1)
#println(translate_8p_ind_to_vec_5(i2,U1,U2,k1,U3,U4,X2,F) - test1)


function translate_8p_vec_to_ind_bSVD(index::Int)

    #u2,K2,U1,s1,K1,U2,E1,X1,s3,K3,U3,E3,X2,U4,E

    temp = index

    x2 = 0
    U1 = 0
    i1 = 0
    U2 = 0
    X1 = 0
    i3 = 0
    U3 = 0
    X2 = 0
    F = 0

    while temp > 0

        x2 += 1

        temp -= len_i_U_i_U_X_i_U_X_F[x2]

    end

    temp += len_i_U_i_U_X_i_U_X_F[x2]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_i_U_X_F[x2,x2,U1]

    end

    temp += len_U_i_U_X_i_U_X_F[x2,x2,U1]

    (a1,b1) = SU[x2,x2,U1,2:3]

    while temp > 0

        i1 += 1

        temp -= len_i_U_X_i_U_X_F[a1,b1,i1]

    end

    temp += len_i_U_X_i_U_X_F[a1,b1,i1]

    while temp > 0

        U2 += 1

        temp -= len_U_X_i_U_X_F[a1,b1,i1,i1,U2]

    end

    temp += len_U_X_i_U_X_F[a1,b1,i1,i1,U2]

    (x1,y1) = SU[i1,i1,U2,2:3]

    while temp > 0

        X1 += 1

        temp -= len_X_i_U_X_F[a1,b1,x1,y1,X1]

    end

    temp += len_X_i_U_X_F[a1,b1,x1,y1,X1]

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    while temp > 0

        i3 += 1

        temp -= len_i_U_X_F[a2,b2,i3]

    end

    temp += len_i_U_X_F[a2,b2,i3]

    while temp > 0

        U3 += 1

        temp -= len_U_X_F[a2,b2,i3,i3,U3]

    end

    temp += len_U_X_F[a2,b2,i3,i3,U3]

    (x3,y3) = SU[i3,i3,U3,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_F[a2,b2,x3,y3,X2]

    end

    temp += len_X_F[a2,b2,x3,y3,X2]

    F = temp

    return (x2,U1,i1,U2,X1,i3,U3,X2,F)

end


function translate_8p_vec_to_ind_naive_bSVD(index::Int)

    count = 0

    x20 = 0
    U10 = 0
    i10 = 0
    U20 = 0
    X10 = 0
    i30 = 0
    U30 = 0
    X20 = 0
    F0 = 0

    for x2 in 1:k+1

        for U1 in 1:RangeU[x2,x2]

            (i5,M1,N1) = SU[x2,x2,U1,:]

            for i1 in 1:k+1

                for U2 in 1:RangeU[i1,i1]

                    (k1,x1,y1) = SU[i1,i1,U2,:]

                    for X1 in 1:RangeX[M1,N1,x1,y1]

                        (M5,N5) = SX[M1,N1,x1,y1,X1,:]

                        for i3 in 1:k+1

                            for U3 in 1:RangeU[i3,i3]

                                (k3,x3,y3) = SU[i3,i3,U3,:]

                                for X2 in 1:RangeX[M5,N5,x3,y3]

                                    (M3,N3) = SX[M5,N5,x3,y3,X2,:]

                                    for F in 1:RangeF[M3,N3]

                                        (k6,x4) = SF[M3,N3,F,:]

                                        count += 1

                                        if count == index

                                            (x20,U10,i10,U20,X10,i30,U30,X20,F0) =
                                                (x2,U1,i1,U2,X1,i3,U3,X2,F)

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

    return (x20,U10,i10,U20,X10,i30,U30,X20,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_8p_ind_to_vec_bSVD(x2,U1,i1,U2,X1,i3,U3,X2,F)

    temp = 0

    if x2 > 1

        temp += sum(len_i_U_i_U_X_i_U_X_F[1:x2-1])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_i_U_X_F[x2,x2,1:U1-1])

    end

    (a1,b1) = SU[x2,x2,U1,2:3]

    if i1 > 1

        temp += sum(len_i_U_X_i_U_X_F[a1,b1,1:i1-1])

    end

    if U2 > 1

        temp += sum(len_U_X_i_U_X_F[a1,b1,i1,i1,1:U2-1])

    end

    (x1,y1) = SU[i1,i1,U2,2:3]

    if X1 > 1

        temp += sum(len_X_i_U_X_F[a1,b1,x1,y1,1:X1-1])

    end

    (a2,b2) = SX[a1,b1,x1,y1,X1,:]

    if i3 > 1

        temp += sum(len_i_U_X_F[a2,b2,1:i3-1])

    end

    if U3 > 1

        temp += sum(len_U_X_F[a2,b2,i3,i3,1:U3-1])

    end

    (x3,y3) = SU[i3,i3,U3,2:3]

    if X2 > 1

        temp += sum(len_X_F[a2,b2,x3,y3,1:X2-1])

    end

    temp += F

    return temp

end


#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(x21,U11,i11,U21,X11,i31,U31,X21,F1) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(x21,U11,i11,U21,X11,i31,U31,X21,F1) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(x21,U11,i11,U21,X11,i31,U31,X21,F1) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(x21,U11,i11,U21,X11,i31,U31,X21,F1) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(x21,U11,i11,U21,X11,i31,U31,X21,F1) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(x21,U11,i11,U21,X11,i31,U31,X21,F1) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(x21,U11,i11,U21,X11,i31,U31,X21,F1) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(x21,U11,i11,U21,X11,i31,U31,X21,F1) - test1)

#test1 = rand(1:configs_8p)
#println(@time translate_8p_vec_to_ind_bSVD(test1))
#println(@time translate_8p_vec_to_ind_naive_bSVD(test1))
#(x21,U11,i11,U21,X11,i31,U31,X21,F1) = translate_8p_vec_to_ind_bSVD(test1)
#println(@time translate_8p_ind_to_vec_bSVD(x21,U11,i11,U21,X11,i31,U31,X21,F1) - test1)
