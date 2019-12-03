function glueing_6p_states_alt_vec(state_6p_1::Array{Complex{Float64},1},state_6p_2::Array{Complex{Float64},1})

    newtensor = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,U5,U6,U7,U8) = translate_10p_vec_to_ind_1(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        (k6,l6,c3,d3) = SU[i6,j6,U5,:]

        (k4,l4,c2,d2) = SU[c3,d3,U6,:]

        (k2,l2,c1,d1) = SU[c2,d2,U7,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        T1 = Uf[k1,l1,k3,l3,c1,d1]

        T2 = Uf[c1,d1,k2,l2,c2,d2]

        T3 = Uf[c2,d2,k4,l4,c3,d3]

        T4 = Uf[c3,d3,k6,l6,i6,j6]

        index_left_old = translate_6p_ind_to_vec(i1,j1,U1,U2,U3,U4)
        index_right_old = translate_6p_ind_to_vec(k1,l1,T1,T2,T3,T4)

        newtensor[index] =
            state_6p_1[index_left_old] *
            state_6p_2[index_right_old] / v(i6) / v(j6) * D^2

    end

    return chop.(newtensor)

end


function ten_p_trafo_1_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,U5,U6,U7,U8) = translate_10p_vec_to_ind_1(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (k6,l6,c4,d4) = SU[M1,N1,U4,:]

        (i4,j4,c3,d3) = SU[c4,d4,U5,:]

        (k4,l4,c2,d2) = SU[c3,d3,U6,:]

        (k2,l2,c1,d1) = SU[c2,d2,U7,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T4 in 1:RangeU[M1,N1]

            (i4old,j4old,i6,j6) = SU[M1,N1,T4,:]

            if  i4 == i4old &&
                j4 == j4old &&
                coupling_rules(i6,k6,c3) == 1 &&
                coupling_rules(j6,l6,d3) == 1

                T5 = Uf[i6,j6,k6,l6,c3,d3]

                index_old = translate_10p_ind_to_vec_1(i1,j1,U1,U2,U3,T4,T5,U6,U7,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i4,M1,i6,k6,c3,c4) * sixjr(j4,N1,j6,l6,d3,d4)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_2_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,i4,j4,Z4,X4,U7,U8) = translate_10p_vec_to_ind_2(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (k6,l6,c4,d4) = SU[M1,N1,U4,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (c2,d2) = SX[c4,d4,x4,y4,X4,:]

        (k2,l2,c1,d1) = SU[c2,d2,U7,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T5 in 1:RangeU[c4,d4]

            (i4old,j4old,c3,d3) = SU[c4,d4,T5,:]

            if  i4 == i4old &&
                j4 == j4old &&
                coupling_rules(c3,k4,c2) == 1 &&
                coupling_rules(d3,l4,d2) == 1

                T6 = Uf[c3,d3,k4,l4,c2,d2]

                index_old = translate_10p_ind_to_vec_1(i1,j1,U1,U2,U3,U4,T5,T6,U7,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i4,c4,c3,c2,k4,x4) * sixjr(j4,d4,d3,d2,l4,y4)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_3_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,i4,j4,Z4,U7,X4,U8) = translate_10p_vec_to_ind_3(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (k6,l6,c4,d4) = SU[M1,N1,U4,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (k2,l2,m2,n2) = SU[x4,y4,U7,:]

        (c1,d1) = SX[c4,d4,m2,n2,X4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for X4old in 1:RangeX[c4,d4,x4,y4]

            (c2,d2) = SX[c4,d4,x4,y4,X4old,:]

            if coupling_rules(c2,k2,c1) == 1 &&
                coupling_rules(d2,l2,d1) == 1

                T7 = Uf[c2,d2,k2,l2,c1,d1]

                index_old = translate_10p_ind_to_vec_2(i1,j1,U1,U2,U3,U4,i4,j4,Z4,X4old,T7,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(x4,c4,c2,c1,k2,m2) * sixjr(y4,d4,d2,d1,l2,n2) *
                    conj(R_matrix(k2,x4,m2)) * R_matrix(l2,y4,n2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_4_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,U5,i4,j4,Z4,X4,U8) =
        translate_10p_vec_to_ind_4(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (k6,l6,c4,d4) = SU[M1,N1,U4,:]

        (k2,l2,c2,d2) = SU[c4,d4,U5,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T7 in 1:RangeU[x4,y4]

            (k2old,l2old,m2,n2) = SU[x4,y4,T7,:]

            if  k2 == k2old &&
                l2 == l2old &&
                coupling_rules(c4,m2,c1) == 1 &&
                coupling_rules(d4,n2,d1) == 1

                X4old = Xf[c4,d4,m2,n2,c1,d1]

                index_old = translate_10p_ind_to_vec_3(i1,j1,U1,U2,U3,U4,i4,j4,Z4,T7,X4old,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(x4,k2,m2,c4,c1,c2) * sixjr(y4,l2,n2,d4,d1,d2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_5_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,U5,i4,j4,Z4,X4,U8) =
        translate_10p_vec_to_ind_4(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (k2,l2,M2,N2) = SU[M1,N1,U4,:]

        (k6,l6,c2,d2) = SU[M2,N2,U5,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T4 in 1:RangeU[M1,N1]

            (k6old,l6old,c4,d4) = SU[M1,N1,T4,:]

            if  k6 == k6old &&
                l6 == l6old &&
                coupling_rules(c4,k2,c2) == 1 &&
                coupling_rules(d4,l2,d2) == 1

                T5 = Uf[c4,d4,k2,l2,c2,d2]

                index_old = translate_10p_ind_to_vec_4(i1,j1,U1,U2,U3,T4,T5,i4,j4,Z4,X4,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(M1,k6,c4,c2,k2,M2) * sixjr(N1,l6,d4,d2,l2,N2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_6_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,j1,U1,U2,U3,U4,U5,i4,j4,Z4,X4,U8) =
        translate_10p_vec_to_ind_4(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (k2,l2,a3,b3) = SU[a2,b2,U3,:]

        (i5,j5,M2,N2) = SU[a3,b3,U4,:]

        (k6,l6,c2,d2) = SU[M2,N2,U5,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T3 in 1:RangeU[a2,b2]

            (i5old,j5old,M1,N1) = SU[a2,b2,T3,:]

            if  i5 == i5old &&
                j5 == j5old &&
                coupling_rules(M1,M2,k2) == 1 &&
                coupling_rules(N1,N2,l2) == 1

                T4 = Uf[M1,N1,k2,l2,M2,N2]

                index_old = translate_10p_ind_to_vec_4(i1,j1,U1,U2,T3,T4,U5,i4,j4,Z4,X4,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(a2,i5,M1,M2,k2,a3) * sixjr(b2,j5,N1,N2,l2,b3)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_bSVD_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p_bSVD
    for index in 1:configs_10p

        (i1,j1,U1,i2,j2,Z2,X2,U4,U5,i4,j4,Z4,X4,U8) = translate_10p_vec_to_ind_bSVD(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (k2,l2,x2,y2) = SU[i2,j2,Z2,:]

        (a3,b3) = SX[a1,b1,x2,y2,X2,:]

        (i5,j5,M2,N2) = SU[a3,b3,U4,:]

        (k6,l6,c2,d2) = SU[M2,N2,U5,:]

        (k4,l4,x4,y4) = SU[i4,j4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U8,:]

        for T2 in 1:RangeU[a1,b1]

            (i2old,j2old,a2,b2) = SU[a1,b1,T2,:]

            if i2 == i2old &&
                j2 == j2old &&
                coupling_rules(a2,k2,a3) == 1 &&
                coupling_rules(b2,l2,b3) == 1

                T3 = Uf[a2,b2,k2,l2,a3,b3]

                index_old = translate_10p_ind_to_vec_4(i1,j1,U1,T2,T3,U4,U5,i4,j4,Z4,X4,U8)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i2,a1,a2,a3,k2,x2) * sixjr(j2,b1,b2,b3,l2,y2) /
                    v(x2) / v(y2) / v(x4) / v(y4) / D^2

            end

        end

    end

    return chop.(Amp_new)

end
