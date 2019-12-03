function transform_6p_state_end_1_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,j5,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i1,j1,a1,b1) = SU[i5,j5,U1,:]

        (i2,j2,M5,N5) = SU[a1,b1,U2,:]

        (i3,j3,M3,N3) = SU[M5,N5,U3,:]

        (i6,j6,i4,j4) = SU[M3,N3,U4,:]

        for T1 in 1:RangeU[i2,j2]

            (i5old,j5old,M1,N1) = SU[i2,j2,T1,:]

            if  i5 == i5old &&
                j5 == j5old &&
                coupling_rules(M1,i1,M5) == 1 &&
                coupling_rules(N1,j1,N5) == 1

                T2 = Uf[M1,N1,i1,j1,M5,N5]

                index_old = translate_6p_ind_to_vec(i2,j2,T1,T2,U3,U4)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(i5,i2,M1,M5,i1,a1) * sixjr(j5,j2,N1,N5,j1,b1) #*
                                    #conj(R_matrix(i4,j4,s4)) * R_matrix(i6,j6,s6)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_end_2_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,j5,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i1,j1,a1,b1) = SU[i5,j5,U1,:]

        (i2,j2,M5,N5) = SU[a1,b1,U2,:]

        (i4,j4,a4,b4) = SU[M5,N5,U3,:]

        (i6,j6,i3,j3) = SU[a4,b4,U4,:]

        for T3 in 1:RangeU[M5,N5]

            (i3old,j3old,M3,N3) = SU[M5,N5,T3,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(M3,i6,i4) == 1 &&
                coupling_rules(N3,j6,j4) == 1

                T4 = Uf[M3,N3,i6,j6,i4,j4]

                index_old = translate_6p_ind_to_vec(i5,j5,U1,U2,T3,T4)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(i3,M5,M3,i4,i6,a4) * sixjr(j3,N5,N3,j4,j6,b4) *
                    conj(R_matrix(i6,i3,a4)) * R_matrix(j6,j3,b4)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_end_3_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,j5,U1,i6,j6,U2,X,U3) = translate_6p_end_vec_to_ind(index)

        (i1,j1,a1,b1) = SU[i5,j5,U1,:]

        (i3,j3,a4,b4) = SU[i6,j6,U2,:]

        (a3,b3) = SX[a1,b1,a4,b4,X,:]

        (i2,j2,i4,j4) = SU[a3,b3,U3,:]

        for T2 in 1:RangeU[a1,b1]

            (i2old,j2old,M5,N5) = SU[a1,b1,T2,:]

            if  i2 == i2old &&
                j2 == j2old &&
                coupling_rules(M5,a4,i4) == 1 &&
                coupling_rules(N5,b4,j4) == 1

                T3 = Uf[M5,N5,i4,j4,a4,b4]
                T4 = Uf[a4,b4,i6,j6,i3,j3]

                index_old = translate_6p_ind_to_vec(i5,j5,U1,T2,T3,T4)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(a1,i2,M5,i4,a4,a3) * sixjr(b1,j2,N5,j4,b4,b3)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_end_4_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,j5,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i1,j1,a1,b1) = SU[i5,j5,U1,:]

        (i6,j6,a2,b2) = SU[a1,b1,U2,:]

        (i3,j3,a3,b3) = SU[a2,b2,U3,:]

        (i2,j2,i4,j4) = SU[a3,b3,U4,:]

        for T2 in 1:RangeU[i6,j6]

            (i3old,j3old,a4,b4) = SU[i6,j6,T2,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(a1,a4,a3) == 1 &&
                coupling_rules(b1,b4,b3) == 1

                X = Xf[a1,b1,a4,b4,a3,b3]

                index_old = translate_6p_end_ind_to_vec(i5,j5,U1,i6,j6,T2,X,U4)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(a1,a3,a4,i3,i6,a2) * sixjr(b1,b3,b4,j3,j6,b2)

            end

        end

    end

    return chop.(Amp_new)

end
