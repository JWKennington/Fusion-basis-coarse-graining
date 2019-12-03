function transform_6p_state_1_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i5,i6) = SF[a3,b3,F,:]

        for T3 in 1:RangeU[a2,b2]

            (i5old,M1,N1) = SU[a2,b2,T3,:]

            if i5 == i5old &&
                coupling_rules(M1,i4,i6) == 1 &&
                coupling_rules(N1,i4,i6) == 1

                Ftild = Ff[M1,N1,i4,i6]

                index_old = translate_6p_ind_to_vec(i1,U1,U2,T3,Ftild)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i4,a2,a3,i5,i6,M1) * sixjr(i4,b2,b3,i5,i6,N1)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_2_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,i6) = SF[M1,N1,F,:]

        for T1 in 1:RangeU[i2,i2]

            (i3old,m3,n3) = SU[i2,i2,T1,:]

            if i3 == i3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(i1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,T1,T2,U3,F)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i2,i1,a1,a2,i3,m3) * sixjr(i2,i1,b1,b2,i3,n3) *
                    conj(R_matrix(i2,i3,m3)) * R_matrix(i2,i3,n3)

            end

        end

    end

    return chop.(Amp_new)

end


function inverse_tranform_6p_3_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p

        (i2,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,m3,n3) = SU[i2,i2,U1,:]

        (i1,a2,b2) = SU[m3,n3,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,Ki6) = SF[M1,N1,F,:]

        for T1 in 1:RangeU[i1,i1]

            (i3old,a1,b1) = SU[i1,i1,T1,:]

            if  s3 == s3old &&
                coupling_rules(a1,i2,a2) == 1 &&
                coupling_rules(b1,i2,b2) == 1

                T2 = Uf[a1,b1,i2,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,T1,T2,U3,F)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i3,i1,a1,a2,i2,m3) * sixjr(i3,i1,b1,b2,i2,n3)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_other_1_inv_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i6,i5) = SF[a3,b3,F,:]

        for T1 in 1:RangeU[i2,i2]

            (i3old,m3,n3) = SU[i2,i2,T1,:]

            if  i3 == i3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(i1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,T1,T2,U3,F)

                Ampl_new[index] += state_6p_other[index_old] *
                    sixjr(i3,a2,a1,i1,i2,m3) * sixjr(i3,b2,b1,i1,i2,n3) *
                    R_matrix(i2,i3,m3) * conj(R_matrix(i2,i3,n3))

            end

        end

    end

    return chop.(Ampl_new)

end


function transform_6p_state_other_2_inv_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    Threads.@threads for index in 1:configs_6p

        (i2,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,m3,n3) = SU[i2,i2,U1,:]

        (i1,a2,b2) = SU[m3,n3,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i6,i5) = SF[a3,b3,F,:]

        for T1 in 1:RangeU[i1,i1]

            (i3old,a1,b1) = SU[i1,i1,T1,:]

            if  i3 == i3old &&
                coupling_rules(i2,a1,a2) == 1 &&
                coupling_rules(i2,b1,b2) == 1

                T2 = Uf[a1,b1,i2,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,T1,T2,U3,F)

                Ampl_new[index] += state_6p_other[index_old] *
                    sixjr(i3,i2,m3,a2,i1,a1) * sixjr(i3,i2,n3,b2,i1,b1)

            end

        end

    end

    return chop.(Ampl_new)

end
