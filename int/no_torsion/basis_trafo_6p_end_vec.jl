function transform_6p_state_end_1_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i1,a1,b1) = SU[i5,i5,U1,:]

        (i2,M5,N5) = SU[a1,b1,U2,:]

        (i3,M3,N3) = SU[M5,N5,U3,:]

        (i6,i4) = SU[M3,N3,F,:]

        for T1 in 1:RangeU[i2,i2]

            (i5old,M1,N1) = SU[i2,i2,T1,:]

            if i5 == i5old &&
                coupling_rules(M1,i1,M5) == 1 &&
                coupling_rules(N1,i1,N5) == 1

                T2 = Uf[M1,N1,i1,M5,N5]

                index_old = translate_6p_ind_to_vec(i2,T1,T2,U3,F)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(i5,i2,M1,M5,i1,a1) * sixjr(i5,i2,N1,N5,i1,b1) #*
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

        (i5,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i1,a1,b1) = SU[i5,i5,U1,:]

        (i2,M5,N5) = SU[a1,b1,U2,:]

        (i4,a4,b4) = SU[M5,N5,U3,:]

        (i6,i3) = SF[a4,b4,F,:]

        for T3 in 1:RangeU[M5,N5]

            (i3old,M3,N3) = SU[M5,N5,T3,:]

            if i3 == i3old &&
                coupling_rules(M3,i6,i4) == 1 &&
                coupling_rules(N3,i6,i4) == 1

                Fold = Ff[M3,N3,i6,i4]

                index_old = translate_6p_ind_to_vec(i5,U1,U2,T3,Fold)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(i3,M5,M3,i4,i6,a4) * sixjr(i3,N5,N3,i4,i6,b4) *
                    conj(R_matrix(i6,i3,a4)) * R_matrix(i6,i3,b4)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_end_3_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,U1,i6,U2,X,F) = translate_6p_end_vec_to_ind(index)

        (i1,a1,b1) = SU[i5,i5,U1,:]

        (i3,a4,b4) = SU[i6,i6,U2,:]

        (a3,b3) = SX[a1,b1,a4,b4,X,:]

        (i2,i4) = SF[a3,b3,F,:]

        for T2 in 1:RangeU[a1,b1]

            (i2old,M5,N5) = SU[a1,b1,T2,:]

            if i2 == i2old &&
                coupling_rules(M5,a4,i4) == 1 &&
                coupling_rules(N5,b4,i4) == 1

                T3 = Uf[M5,N5,i4,a4,b4]
                Fold = Ff[a4,b4,i6,i3]

                index_old = translate_6p_ind_to_vec(i5,U1,T2,T3,Fold)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(a1,i2,M5,i4,a4,a3) * sixjr(b1,i2,N5,i4,b4,b3)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_end_4_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i5,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i1,a1,b1) = SU[i5,i5,U1,:]

        (i6,a2,b2) = SU[a1,b1,U2,:]

        (i3,a3,b3) = SU[a2,b2,U3,:]

        (i2,i4) = SU[a3,b3,F,:]

        for T2 in 1:RangeU[i6,i6]

            (i3old,a4,b4) = SU[i6,i6,T2,:]

            if i3 == i3old &&
                coupling_rules(a1,a4,a3) == 1 &&
                coupling_rules(b1,b4,b3) == 1

                X = Xf[a1,b1,a4,b4,a3,b3]

                index_old = translate_6p_end_ind_to_vec(i5,U1,i6,T2,X,F)

                Amp_new[index] +=
                    state_6p[index_old] *
                    sixjr(a1,a3,a4,i3,i6,a2) * sixjr(b1,b3,b4,i3,i6,b2)

            end

        end

    end

    return chop.(Amp_new)

end
