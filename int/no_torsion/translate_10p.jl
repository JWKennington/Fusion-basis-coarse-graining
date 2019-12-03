#Define the functions to translate back and forth between 6p and 10p states

function translate_10p_vec_to_ind_1(index::Int)

    temp = index

    i = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    U5 = 0
    U6 = 0
    U7 = 0
    F = 0

    while temp > 0

        i += 1

        temp -= len_i_7U_F[i]

    end

    temp += len_i_7U_F[i]

    while temp > 0

        U1 += 1

        temp -= len_7U_F[i,i,U1]

    end

    temp += len_7U_F[i,i,U1]

    (a1,b1) = SU[i,i,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_6U_F[a1,b1,U2]

    end

    temp += len_6U_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_5U_F[a2,b2,U3]

    end

    temp += len_5U_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_4U_F[a3,b3,U4]

    end

    temp += len_4U_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    while temp > 0

        U5 += 1

        temp -= len_3U_F[a4,b4,U5]

    end

    temp += len_3U_F[a4,b4,U5]

    (a5,b5) = SU[a4,b4,U5,2:3]

    while temp > 0

        U6 += 1

        temp -= len_2U_F[a5,b5,U6]

    end

    temp += len_2U_F[a5,b5,U6]

    (a6,b6) = SU[a5,b5,U6,2:3]

    while temp > 0

        U7 += 1

        temp -= len_U_F[a6,b6,U7]

    end

    temp += len_U_F[a6,b6,U7]

    F = temp

    return (i,U1,U2,U3,U4,U5,U6,U7,F)

end


function translate_10p_vec_to_ind_naive_1(index::Int)

    count = 0

    i0 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    U50 = 0
    U60 = 0
    U70 = 0
    F0 = 0

    for i in 1:k+1

        for U1 in 1:RangeU[i,i]

            (i3,a1,b1) = SU[i,i,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i2,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i5,M1,N1) = SU[a2,b2,U3,:]

                    for U4 in 1:RangeU[M1,N1]

                        (i4,i6,j6) = SU[M1,N1,U4,:]

                        for U5 in 1:RangeU[i6,j6]

                            (k6,c3,d3) = SU[i6,j6,U5,:]

                            for U6 in 1:RangeU[c3,d3]

                                (k4,c2,d2) = SU[c3,d3,U6,:]

                                for U7 in 1:RangeU[c2,d2]

                                    (k2,c1,d1) = SU[c2,d2,U7,:]

                                    for F in 1:RangeF[c1,d1]

                                        count += 1

                                        if count == index

                                            (i0,U10,U20,U30,U40,U50,U60,U70,F0) =
                                                (i,U1,U2,U3,U4,U5,U6,U7,F)

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

    return (i0,U10,U20,U30,U40,U50,U60,U70,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_10p_ind_to_vec_1(i::Int,U1::Int,U2::Int,U3::Int,
    U4::Int,U5::Int,U6::Int,U7::Int,F::Int)

    temp = 0

    if i > 1

        temp += sum(len_i_7U_F[1:i-1])

    end

    if U1 > 1

        temp += sum(len_7U_F[i,i,1:U1-1])

    end

    (a1,b1) = SU[i,i,U1,2:3]

    if U2 > 1

        temp += sum(len_6U_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_5U_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_4U_F[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if U5 > 1

        temp += sum(len_3U_F[a4,b4,1:U5-1])

    end

    (a5,b5) = SU[a4,b4,U5,2:3]

    if U6 > 1

        temp += sum(len_2U_F[a5,b5,1:U6-1])

    end

    (a6,b6) = SU[a5,b5,U6,2:3]

    if U7 > 1

        temp += sum(len_U_F[a6,b6,1:U7-1])

    end

    temp += F

    return temp

end




function translate_10p_vec_to_ind_2(index::Int)

    temp = index

    i1 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    i4 = 0
    Z4 = 0
    X4 = 0
    U7 = 0
    F = 0

    while temp > 0

        i1 += 1

        temp -= len_i_4U_i_U_X_U_F[i1]

    end

    temp += len_i_4U_i_U_X_U_F[i1]

    while temp > 0

        U1 += 1

        temp -= len_4U_i_U_X_U_F[i1,i1,U1]

    end

    temp += len_4U_i_U_X_U_F[i1,i1,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_3U_i_U_X_U_F[a1,b1,U2]

    end

    temp += len_3U_i_U_X_U_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_2U_i_U_X_U_F[a2,b2,U3]

    end

    temp += len_2U_i_U_X_U_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_U_i_U_X_U_F[a3,b3,U4]

    end

    temp += len_U_i_U_X_U_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    while temp > 0

        i4 += 1

        temp -= len_i_U_X_U_F[a4,b4,i4]

    end

    temp += len_i_U_X_U_F[a4,b4,i4]

    while temp > 0

        Z4 += 1

        temp -= len_U_X_U_F[a4,b4,i4,i4,Z4]

    end

    temp += len_U_X_U_F[a4,b4,i4,i4,Z4]

    (a7,b7) = SU[i4,i4,Z4,2:3]

    while temp > 0

        X4 += 1

        temp -= len_X_U_F[a4,b4,a7,b7,X4]

    end

    temp += len_X_U_F[a4,b4,a7,b7,X4]

    (a8,b8) = SX[a4,b4,a7,b7,X4,:]

    while temp > 0

        U7 += 1

        temp -= len_U_F[a8,b8,U7]

    end

    temp += len_U_F[a8,b8,U7]

    F = temp

    return (i1,U1,U2,U3,U4,i4,Z4,X4,U7,F)

end


function translate_10p_vec_to_ind_naive_2(index::Int)

    count = 0

    i10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    i40 = 0
    Z40 = 0
    X40 = 0
    U70 = 0
    F0 = 0

    for i1 in 1:k+1

        for U1 in 1:RangeU[i1,i1]

            (i3,a1,b1) = SU[i1,i1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i2,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i5,M1,N1) = SU[a2,b2,U3,:]

                    for U4 in 1:RangeU[M1,N1]

                        (k6,c4,d4) = SU[M1,N1,U4,:]

                        for i4 in 1:k+1

                            for Z4 in 1:RangeU[i4,i4]

                                (k4,x4,y4) = SU[i4,i4,Z4,:]

                                for X4 in 1:RangeX[c4,d4,x4,y4]

                                    (c2,d2) = SX[c4,d4,x4,y4,X4,:]

                                    for U7 in 1:RangeU[c2,d2]

                                        (k2,c1,d1) = SU[c2,d2,U7,:]

                                        for F in 1:RangeF[c1,d1]

                                            count += 1

                                            if count == index

                                                (i10,U10,U20,U30,U40,i40,Z40,X40,U70,F0) =
                                                    (i1,U1,U2,U3,U4,i4,Z4,X4,U7,F)

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

    end

    return (i10,U10,U20,U30,U40,i40,Z40,X40,U70,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_10p_ind_to_vec_2(i1::Int,U1::Int,U2::Int,
    U3::Int,U4::Int,
    i4::Int,Z4::Int,X4::Int,U7::Int,F::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_4U_i_U_X_U_F[1:i1-1])

    end

    if U1 > 1

        temp += sum(len_4U_i_U_X_U_F[i1,i1,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if U2 > 1

        temp += sum(len_3U_i_U_X_U_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_2U_i_U_X_U_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_U_i_U_X_U_F[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if i4 > 1

        temp += sum(len_i_U_X_U_F[a4,b4,1:i4-1])

    end

    if Z4 > 1

        temp += sum(len_U_X_U_F[a4,b4,i4,i4,1:Z4-1])

    end

    (a5,b5) = SU[i4,i4,Z4,2:3]

    if X4 > 1

        temp += sum(len_X_U_F[a4,b4,a5,b5,1:X4-1])

    end

    (a6,b6) = SX[a4,b4,a5,b5,X4,:]

    if U7 > 1

        temp += sum(len_U_F[a6,b6,1:U7-1])

    end

    temp += F

    return temp

end



# 3rd translation between 10p indices and vectorized form


function translate_10p_vec_to_ind_3(index::Int)

    temp = index

    i1 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    i4 = 0
    Z4 = 0
    U7 = 0
    X4 = 0
    F = 0

    while temp > 0

        i1 += 1

        temp -= len_i_4U_i_2U_X_F[i1]

    end

    temp += len_i_4U_i_2U_X_F[i1]

    while temp > 0

        U1 += 1

        temp -= len_4U_i_2U_X_F[i1,i1,U1]

    end

    temp += len_4U_i_2U_X_F[i1,i1,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_3U_i_2U_X_F[a1,b1,U2]

    end

    temp += len_3U_i_2U_X_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_2U_i_2U_X_F[a2,b2,U3]

    end

    temp += len_2U_i_2U_X_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_U_i_2U_X_F[a3,b3,U4]

    end

    temp += len_U_i_2U_X_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    while temp > 0

        i4 += 1

        temp -= len_i_2U_X_F[a4,b4,i4]

    end

    temp += len_i_2U_X_F[a4,b4,i4]

    while temp > 0

        Z4 += 1

        temp -= len_2U_X_F[a4,b4,i4,i4,Z4]

    end

    temp += len_2U_X_F[a4,b4,i4,i4,Z4]

    (a7,b7) = SU[i4,i4,Z4,2:3]

    while temp > 0

        U7 += 1

        temp -= len_U_X_F[a4,b4,a7,b7,U7]

    end

    temp += len_U_X_F[a4,b4,a7,b7,U7]

    (a8,b8) = SU[a7,b7,U7,2:3]

    while temp > 0

        X4 += 1

        temp -= len_X_F[a4,b4,a8,b8,X4]

    end

    temp += len_X_F[a4,b4,a8,b8,X4]

    F = temp

    return (i1,U1,U2,U3,U4,i4,Z4,U7,X4,F)

end


function translate_10p_vec_to_ind_naive_3(index::Int)

    count = 0

    i10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    i40 = 0
    Z40 = 0
    U70 = 0
    X40 = 0
    F0 = 0

    for i1 in 1:k+1

        for U1 in 1:RangeU[i1,i1]

            (i3,a1,b1) = SU[i1,i1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i2,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i5,M1,N1) = SU[a2,b2,U3,:]

                    for U4 in 1:RangeU[M1,N1]

                        (k6,c4,d4) = SU[M1,N1,U4,:]

                        for i4 in 1:k+1

                            for Z4 in 1:RangeU[i4,i4]

                                (k4,x4,y4) = SU[i4,i4,Z4,:]

                                for U7 in 1:RangeU[x4,y4]

                                    (k2,m2,n2) = SU[x4,y4,U7,:]

                                    for X4 in 1:RangeX[c4,d4,m2,n2]

                                        (c1,d1) = SX[c4,d4,m2,n2,X4,:]

                                        for F in 1:RangeF[c1,d1]

                                            count += 1

                                            if count == index

                                                (i10,U10,U20,U30,U40,i40,Z40,U70,X40,F0) =
                                                    (i1,U1,U2,U3,U4,i4,Z4,U7,X4,F)

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

    end

    return (i10,U10,U20,U30,U40,i40,Z40,U70,X40,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_10p_ind_to_vec_3(i1::Int,U1::Int,U2::Int,
    U3::Int,U4::Int,
    i4::Int,Z4::Int,U7::Int,X4::Int,F::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_4U_i_2U_X_F[1:i1-1])

    end

    if U1 > 1

        temp += sum(len_4U_i_2U_X_F[i1,i1,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if U2 > 1

        temp += sum(len_3U_i_2U_X_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_2U_i_2U_X_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_U_i_2U_X_F[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if i4 > 1

        temp += sum(len_i_2U_X_F[a4,b4,1:i4-1])

    end

    if Z4 > 1

        temp += sum(len_2U_X_F[a4,b4,i4,i4,1:Z4-1])

    end

    (a5,b5) = SU[i4,i4,Z4,2:3]

    if U7 > 1

        temp += sum(len_U_X_F[a4,b4,a5,b5,1:U7-1])

    end

    (a6,b6) = SU[a5,b5,U7,2:3]

    if X4 > 1

        temp += sum(len_X_F[a4,b4,a6,b6,1:X4-1])

    end

    (a7,b7) = SX[a4,b4,a6,b6,X4,:]

    temp += F

    return temp

end



# 4th translation between 10p indices and vectorized form


function translate_10p_vec_to_ind_4(index::Int)

    temp = index

    i1 = 0
    U1 = 0
    U2 = 0
    U3 = 0
    U4 = 0
    U5 = 0
    i4 = 0
    Z4 = 0
    X4 = 0
    F = 0

    while temp > 0

        i1 += 1

        temp -= len_i_5U_i_U_X_F[i1]

    end

    temp += len_i_5U_i_U_X_F[i1]

    while temp > 0

        U1 += 1

        temp -= len_5U_i_U_X_F[i1,i1,U1]

    end

    temp += len_5U_i_U_X_F[i1,i1,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        U2 += 1

        temp -= len_4U_i_U_X_F[a1,b1,U2]

    end

    temp += len_4U_i_U_X_F[a1,b1,U2]

    (a2,b2) = SU[a1,b1,U2,2:3]

    while temp > 0

        U3 += 1

        temp -= len_3U_i_U_X_F[a2,b2,U3]

    end

    temp += len_3U_i_U_X_F[a2,b2,U3]

    (a3,b3) = SU[a2,b2,U3,2:3]

    while temp > 0

        U4 += 1

        temp -= len_2U_i_U_X_F[a3,b3,U4]

    end

    temp += len_2U_i_U_X_F[a3,b3,U4]

    (a4,b4) = SU[a3,b3,U4,2:3]

    while temp > 0

        U5 += 1

        temp -= len_U_i_U_X_F[a4,b4,U5]

    end

    temp += len_U_i_U_X_F[a4,b4,U5]

    (a5,b5) = SU[a4,b4,U5,2:3]

    while temp > 0

        i4 += 1

        temp -= len_i_U_X_F[a5,b5,i4]

    end

    temp += len_i_U_X_F[a5,b5,i4]

    while temp > 0

        Z4 += 1

        temp -= len_U_X_F[a5,b5,i4,i4,Z4]

    end

    temp += len_U_X_F[a5,b5,i4,i4,Z4]

    (a7,b7) = SU[i4,i4,Z4,2:3]

    while temp > 0

        X4 += 1

        temp -= len_X_F[a5,b5,a7,b7,X4]

    end

    temp += len_X_F[a5,b5,a7,b7,X4]

    F = temp

    return (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F)

end


function translate_10p_vec_to_ind_naive_4(index::Int)

    count = 0

    i10 = 0
    U10 = 0
    U20 = 0
    U30 = 0
    U40 = 0
    U50 = 0
    i40 = 0
    Z40 = 0
    X40 = 0
    F0 = 0

    for i1 in 1:k+1

        for U1 in 1:RangeU[i1,i1]

            (i3,a1,b1) = SU[i1,i1,U1,:]

            for U2 in 1:RangeU[a1,b1]

                (i2,a2,b2) = SU[a1,b1,U2,:]

                for U3 in 1:RangeU[a2,b2]

                    (i5,M1,N1) = SU[a2,b2,U3,:]

                    for U4 in 1:RangeU[M1,N1]

                        (k6,c4,d4) = SU[M1,N1,U4,:]

                        for U5 in 1:RangeU[c4,d4]

                            (k2,c2,d2) = SU[c4,d4,U5,:]

                            for i4 in 1:k+1

                                for Z4 in 1:RangeU[i4,i4]

                                    (k4,x4,y4) = SU[i4,i4,Z4,:]

                                    for X4 in 1:RangeX[c2,d2,x4,y4]

                                        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                                        for F in 1:RangeF[c1,d1]

                                            count += 1

                                            if count == index

                                                (i10,U10,U20,U30,U40,U50,i40,Z40,X40,F0) =
                                                    (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F)

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

    end

    return (i10,U10,U20,U30,U40,U50,i40,Z40,X40,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_10p_ind_to_vec_4(i1::Int,U1::Int,U2::Int,
    U3::Int,U4::Int,U5::Int,
    i4::Int,Z4::Int,X4::Int,F::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_5U_i_U_X_F[1:i1-1])

    end

    if U1 > 1

        temp += sum(len_5U_i_U_X_F[i1,i1,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if U2 > 1

        temp += sum(len_4U_i_U_X_F[a1,b1,1:U2-1])

    end

    (a2,b2) = SU[a1,b1,U2,2:3]

    if U3 > 1

        temp += sum(len_3U_i_U_X_F[a2,b2,1:U3-1])

    end

    (a3,b3) = SU[a2,b2,U3,2:3]

    if U4 > 1

        temp += sum(len_2U_i_U_X_F[a3,b3,1:U4-1])

    end

    (a4,b4) = SU[a3,b3,U4,2:3]

    if U5 > 1

        temp += sum(len_U_i_U_X_F[a4,b4,1:U5-1])

    end

    (a5,b5) = SU[a4,b4,U5,2:3]

    if i4 > 1

        temp += sum(len_i_U_X_F[a5,b5,1:i4-1])

    end

    if Z4 > 1

        temp += sum(len_U_X_F[a5,b5,i4,i4,1:Z4-1])

    end

    (a6,b6) = SU[i4,i4,Z4,2:3]

    if X4 > 1

        temp += sum(len_X_F[a5,b5,a6,b6,1:X4-1])

    end

    temp += F

    return temp

end


# 4th translation between 10p indices and vectorized form


function translate_10p_vec_to_ind_bSVD(index::Int)

    temp = index

    i1 = 0
    U1 = 0
    i2 = 0
    Z2 = 0
    #E2 = 0
    X2 = 0
    U4 = 0
    U5 = 0
    i4 = 0
    Z4 = 0
    #E4 = 0
    X4 = 0
    F = 0

    while temp > 0

        i1 += 1

        temp -= len_i_U_i_U_X_2U_i_U_X_F[i1]

    end

    temp += len_i_U_i_U_X_2U_i_U_X_F[i1]

    while temp > 0

        U1 += 1

        temp -= len_U_i_U_X_2U_i_U_X_F[i1,i1,U1]

    end

    temp += len_U_i_U_X_2U_i_U_X_F[i1,i1,U1]

    (a1,b1) = SU[i1,i1,U1,2:3]

    while temp > 0

        i2 += 1

        temp -= len_i_U_X_2U_i_U_X_F[a1,b1,i2]

    end

    temp += len_i_U_X_2U_i_U_X_F[a1,b1,i2]

    while temp > 0

        Z2 += 1

        temp -= len_U_X_2U_i_U_X_F[a1,b1,i2,i2,Z2]

    end

    temp += len_U_X_2U_i_U_X_F[a1,b1,i2,i2,Z2]

    (x2,y2) = SU[i2,i2,Z2,2:3]

    while temp > 0

        X2 += 1

        temp -= len_X_2U_i_U_X_F[a1,b1,x2,y2,X2]

    end

    temp += len_X_2U_i_U_X_F[a1,b1,x2,y2,X2]

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    while temp > 0

        U4 += 1

        temp -= len_2U_i_U_X_F[a2,b2,U4]

    end

    temp += len_2U_i_U_X_F[a2,b2,U4]

    (a3,b3) = SU[a2,b2,U4,2:3]

    while temp > 0

        U5 += 1

        temp -= len_U_i_U_X_F[a3,b3,U5]

    end

    temp += len_U_i_U_X_F[a3,b3,U5]

    (a4,b4) = SU[a3,b3,U5,2:3]

    while temp > 0

        i4 += 1

        temp -= len_i_U_X_F[a4,b4,i4]

    end

    temp += len_i_U_X_F[a4,b4,i4]

    while temp > 0

        Z4 += 1

        temp -= len_U_X_F[a4,b4,i4,i4,Z4]

    end

    temp += len_U_X_F[a4,b4,i4,i4,Z4]

    (x4,y4) = SU[i4,i4,Z4,2:3]

    while temp > 0

        X4 += 1

        temp -= len_X_F[a4,b4,x4,y4,X4]

    end

    temp += len_X_F[a4,b4,x4,y4,X4]

    F = temp

    return (i1,U1,i2,Z2,X2,U4,U5,i4,Z4,X4,F)

end


function translate_10p_vec_to_ind_naive_bSVD(index::Int)

    count = 0

    i10 = 0
    U10 = 0
    i20 = 0
    Z20 = 0
    X20 = 0
    U40 = 0
    U50 = 0
    i40 = 0
    Z40 = 0
    X40 = 0
    F0 = 0

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

                                for i4 in 1:k+1

                                    for Z4 in 1:RangeU[i4,i4]

                                        (k4,x4,y4) = SU[i4,i4,Z4,:]

                                        for X4 in 1:RangeX[c2,d2,x4,y4]

                                            (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                                            for F in 1:RangeF[c1,d1]

                                                count += 1

                                                if count == index

                                                    (i10,U10,i20,Z20,X20,U40,U50,i40,Z40,X40,F0) =
                                                        (i1,U1,i2,Z2,X2,U4,U5,i4,Z4,X4,F)

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

        end

    end

    return (i10,U10,i20,Z20,X20,U40,U50,i40,Z40,X40,F0)

end

#Define the translation from set of indices to the vectorized index.

function translate_10p_ind_to_vec_bSVD(i1::Int,U1::Int,
    i2::Int,Z2::Int,X2::Int,
    U4::Int,U5::Int,
    i4::Int,Z4::Int,X4::Int,F::Int)

    temp = 0

    if i1 > 1

        temp += sum(len_i_U_i_U_X_2U_i_U_X_F[1:i1-1])

    end

    if U1 > 1

        temp += sum(len_U_i_U_X_2U_i_U_X_F[i1,i1,1:U1-1])

    end

    (a1,b1) = SU[i1,i1,U1,2:3]

    if i2 > 1

        temp += sum(len_i_U_X_2U_i_U_X_F[a1,b1,1:i2-1])

    end

    if Z2 > 1

        temp += sum(len_U_X_2U_i_U_X_F[a1,b1,i2,i2,1:Z2-1])

    end

    (x2,y2) = SU[i2,i2,Z2,2:3]

    if X2 > 1

        temp += sum(len_X_2U_i_U_X_F[a1,b1,x2,y2,1:X2-1])

    end

    (a2,b2) = SX[a1,b1,x2,y2,X2,:]

    if U4 > 1

        temp += sum(len_2U_i_U_X_F[a2,b2,1:U4-1])

    end

    (a3,b3) = SU[a2,b2,U4,2:3]

    if U5 > 1

        temp += sum(len_U_i_U_X_F[a3,b3,1:U5-1])

    end

    (a4,b4) = SU[a3,b3,U5,2:3]

    if i4 > 1

        temp += sum(len_i_U_X_F[a4,b4,1:i4-1])

    end

    if Z4 > 1

        temp += sum(len_U_X_F[a4,b4,i4,i4,1:Z4-1])

    end

    (x4,y4) = SU[i4,i4,Z4,2:3]

    if X4 > 1

        temp += sum(len_X_F[a4,b4,x4,y4,1:X4-1])

    end

    temp += F

    return temp

end


#test1 = rand(1:configs_10p)
#println(@time translate_10p_vec_to_ind_bSVD(test1))
#println(@time translate_10p_vec_to_ind_naive_bSVD(test1))
#(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) = translate_10p_vec_to_ind_bSVD(test1)
#println(@time translate_10p_ind_to_vec_bSVD(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) - test1)

#test1 = rand(1:configs_10p)
#println(@time translate_10p_vec_to_ind_bSVD(test1))
#println(@time translate_10p_vec_to_ind_naive_bSVD(test1))
#(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) = translate_10p_vec_to_ind_bSVD(test1)
#println(@time translate_10p_ind_to_vec_bSVD(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) - test1)

#test1 = rand(1:configs_10p)
#println(@time translate_10p_vec_to_ind_bSVD(test1))
#println(@time translate_10p_vec_to_ind_naive_bSVD(test1))
#(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) = translate_10p_vec_to_ind_bSVD(test1)
#println(@time translate_10p_ind_to_vec_bSVD(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) - test1)

#test1 = rand(1:configs_10p)
#println(@time translate_10p_vec_to_ind_bSVD(test1))
#println(@time translate_10p_vec_to_ind_naive_bSVD(test1))
#(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) = translate_10p_vec_to_ind_bSVD(test1)
#println(@time translate_10p_ind_to_vec_bSVD(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) - test1)

#test1 = rand(1:configs_10p)
#println(@time translate_10p_vec_to_ind_bSVD(test1))
#println(@time translate_10p_vec_to_ind_naive_bSVD(test1))
#(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) = translate_10p_vec_to_ind_bSVD(test1)
#println(@time translate_10p_ind_to_vec_bSVD(i11,U11,i21,Z21,X21,U41,U51,i41,Z41,X41,F1) - test1)
