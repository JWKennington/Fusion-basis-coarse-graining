function transform_8p_state_1_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,i1,j1,U2,X1,U3,U4,U5) = translate_8p_vec_to_ind_2(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i3,j3,a1,b1) = SU[i1,j1,U2,:]

        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

        (k6,l6,c2,d2) = SU[M2,N2,U3,:]

        (x4,y4,c1,d1) = SU[c2,d2,U4,:]

        (k3,l3,k1,l1) = SU[c1,d1,U5,:]

        for T2 in 1:RangeU[a1,b1]

            (x2old,y2old,a3,b3) = SU[a1,b1,T2,:] #We need to sum over a3,b3

            if  x2 == x2old &&
                y2 == y2old &&
                coupling_rules(a3,i5,M2) == 1 &&
                coupling_rules(b3,j5,N2) == 1

                T3 = Uf[a3,b3,i5,j5,M2,N2]

                index_old = translate_8p_ind_to_vec_1(i1,j1,U2,T2,T3,U3,U4,U5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(a1,x2,a3,i5,M2,M1) * sixjr(b1,y2,b3,j5,N2,N1) *
                    conj(R_matrix(x2,i5,M1)) * R_matrix(y2,j5,N1)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_8p_state_2_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,i1,j1,U2,X1,k1,l1,U3,X2,U4) = translate_8p_vec_to_ind_3(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i3,j3,a1,b1) = SU[i1,j1,U2,:]

        (M2,N2) = SX[M1,N1,a1,b1,X1,:]

        (k3,l3,c1,d1) = SU[k1,l1,U3,:]

        (M3,N3) = SX[M2,N2,c1,d1,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U4,:]

        for T3 in 1:RangeU[M2,N2]

            (k6old,l6old,c2,d2) = SU[M2,N2,T3,:]

            if  k6 == k6old &&
                l6 == l6old &&
                coupling_rules(c2,c1,x4) == 1 &&
                coupling_rules(d2,d1,y4) == 1

                T4 = Uf[c2,d2,x4,y4,c1,d1]

                T5 = Uf[c1,d1,k3,l3,k1,l1]

                index_old = translate_8p_ind_to_vec_2(x2,y2,U1,i1,j1,U2,X1,T3,T4,T5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(M2,k6,c2,x4,c1,M3) * sixjr(N2,l6,d2,y4,d1,N3) *
                    R_matrix(x4,k6,M3) * conj(R_matrix(y4,l6,N3)) #*
                    #conj(R_matrix(x4,y4,u4)) * R_matrix(k6,l6,t6)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_8p_state_3_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,U2,U3,k1,l1,U4,X2,U5) = translate_8p_vec_to_ind_4(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i1,j1,M4,N4) = SU[M1,N1,U2,:]

        (i3,j3,M2,N2) = SU[M4,N4,U3,:]

        (k3,l3,c1,d1) = SU[k1,l1,U4,:]

        (M3,N3) = SX[M2,N2,c1,d1,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U5,:]

        for T2 in 1:RangeU[i1,j1]

            (i3old,j3old,a1,b1) = SU[i1,j1,T2,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(a1,M1,M2) == 1 &&
                coupling_rules(b1,N1,N2) == 1

                X1 = Xf[M1,N1,a1,b1,M2,N2]

                index_old = translate_8p_ind_to_vec_3(x2,y2,U1,i1,j1,T2,X1,k1,l1,U4,X2,U5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(i3,i1,a1,M1,M2,M4) * sixjr(j3,j1,b1,N1,N2,N4)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_8p_state_4_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,U2,k1,l1,U3,U4,X2,U5) = translate_8p_vec_to_ind_5(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i1,j1,M4,N4) = SU[M1,N1,U2,:]

        (k3,l3,c1,d1) = SU[k1,l1,U3,:]

        (i3,j3,c2,d2) = SU[c1,d1,U4,:]

        (M3,N3) = SX[M4,N4,c2,d2,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U5,:]

        for T3 in 1:RangeU[M4,N4]

            (i3old,j3old,M2,N2) = SU[M4,N4,T3,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(M2,c1,M3) == 1 &&
                coupling_rules(N2,d1,N3) == 1

                X2old = Xf[M2,N2,c1,d1,M3,N3]

                index_old = translate_8p_ind_to_vec_4(x2,y2,U1,U2,T3,k1,l1,U3,X2old,U5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(i3,M4,M2,M3,c1,c2) * sixjr(j3,N4,N2,N3,d1,d2)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_8p_state_5_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,U2,i3,j3,U3,U4,X2,U5) = translate_8p_vec_to_ind_5(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i1,j1,M4,N4) = SU[M1,N1,U2,:]

        (k3,l3,x3,y3) = SU[i3,j3,U3,:]

        (k1,l1,c2,d2) = SU[x3,y3,U4,:]

        (M3,N3) = SX[M4,N4,c2,d2,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U5,:]

        for T3 in 1:RangeU[k1,l1]

            (k3old,l3old,c1,d1) = SU[k1,l1,T3,:]

            if  k3 == k3old &&
                l3 == l3old &&
                coupling_rules(c1,i3,c2) == 1 &&
                coupling_rules(d1,j3,d2) == 1

                T4 = Uf[c1,d1,i3,j3,c2,d2]

                index_old = translate_8p_ind_to_vec_5(x2,y2,U1,U2,k1,l1,T3,T4,X2,U5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(i3,c2,c1,k1,k3,x3) * sixjr(j3,d2,d1,l1,l3,y3) *
                    conj(R_matrix(k1,x3,c2)) * R_matrix(l1,y3,d2)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_8p_state_6_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (x2,y2,U1,U2,U3,i3,j3,U4,X2,U5) = translate_8p_vec_to_ind_4(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (i1,j1,M4,N4) = SU[M1,N1,U2,:]

        (k1,l1,M5,N5) = SU[M4,N4,U3,:]

        (k3,l3,x3,y3) = SU[i3,j3,U4,:]

        (M3,N3) = SX[M5,N5,x3,y3,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U5,:]

        for T4 in 1:RangeU[x3,y3]

            (k1old,l1old,c2,d2) = SU[x3,y3,T4,:]

            if  k1 == k1old &&
                l1 == l1old &&
                coupling_rules(M4,c2,M3) == 1&&
                coupling_rules(N4,d2,N3) == 1

                X2old = Xf[M4,N4,c2,d2,M3,N3]

                index_old = translate_8p_ind_to_vec_5(x2,y2,U1,U2,i3,j3,U4,T4,X2old,U5)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(M4,M3,c2,x3,k1,M5) * sixjr(N4,N3,d2,y3,l1,N5)

            end

        end

    end

    return chop.(Amp_new)

end


function state_before_SVD_2_vec(state_8p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p_bSVD
    for index in 1:configs_8p

        (x2,y2,U1,i1,j1,U2,X1,i3,j3,U3,X2,U4) = translate_8p_vec_to_ind_bSVD(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (k1,l1,x1,y1) = SU[i1,j1,U2,:]

        (M5,N5) = SX[M1,N1,x1,y1,X1,:]

        (k3,l3,x3,y3) = SU[i3,j3,U3,:]

        (M3,N3) = SX[M5,N5,x3,y3,X2,:]

        (k6,l6,x4,y4) = SU[M3,N3,U4,:]

        for T2 in 1:RangeU[M1,N1]

            (i1old,j1old,M4,N4) = SU[M1,N1,T2,:]

            if  i1 == i1old &&
                j1 == j1old &&
                coupling_rules(M4,k1,M5) == 1 &
                coupling_rules(N4,l1,N5) == 1

                T3 = Uf[M4,N4,k1,l1,M5,N5]

                index_old = translate_8p_ind_to_vec_4(x2,y2,U1,T2,T3,i3,j3,U3,X2,U4)

                Amp_new[index] +=
                    state_8p[index_old] *
                    sixjr(i1,M1,M4,M5,k1,x1) * sixjr(j1,N1,N4,N5,l1,y1) /
                    v(x1) / v(y1) / v(x3) / v(y3) / D^2

            end

        end

    end

    return chop.(Amp_new)

end
