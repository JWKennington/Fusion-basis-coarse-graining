function define_6p_strong_vec(g::Float64)

    #g = 0.5

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i5,i6) = SF[a3,b3,F,:]

        if a1 == b1 && a2 == b2 && a3 == b3

            Amp_new[index] = v(i1) * v(i2) * v(i3) * v(i4) * v(i5) * v(i6) / D^6 *
                heat_ker(g,i1) * heat_ker(g,i2) * heat_ker(g,i3) * heat_ker(g,i4) * heat_ker(g,i5) * heat_ker(g,i6)

        end

    end

    return chop.(Amp_new/Amp_new[1])

end



function define_6p_state_other_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i6,i5) = SF[a3,b3,F,:]

        Ftild = Ff[a3,b3,i5,i6]

        index_old = translate_6p_ind_to_vec(i1,U1,U2,U3,Ftild)

        Amp_new[index] = state_6p[index_old] *
            conj(R_matrix(i2,i1,a1)) * R_matrix(i2,i1,b1) *
            conj(R_matrix(i3,a1,a2)) * R_matrix(i3,b1,b2) *
            conj(R_matrix(i4,a2,a3)) * R_matrix(i4,b2,b3)

    end

    return chop.(Amp_new)

end



function transform_6p_state_1_vec(state_6p::Array{Complex{Float64},1})

    #How to combine punctures into indices?

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,i6) = SF[M1,N1,F,:]

        for T3 in 1:RangeU[a2,b2]

            (i4old,a3,b3) = SU[a2,b2,T3,:]

            if  i4 == i4old &&
                coupling_rules(a3,i5,i6) == 1 &&
                coupling_rules(b3,i5,i6) == 1

                Fold = Ff[a3,b3,i5,i6]

                index_old = translate_6p_ind_to_vec(i1,U1,U2,T3,Fold)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i4,a2,a3,i5,i6,M1) * sixjr(i4,b2,b3,i5,i6,N1)

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_2_vec(state_6p::Array{Complex{Float64},1})

    #How to combine punctures into indices?

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i2,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,m3,n3) = SU[i2,i2,U1,:]

        (i1,a2,b2) = SU[m3,n3,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,i6) = SF[M1,N1,F,:]

        for T1 in 1:RangeU[i1,i1]

            (i2old,a1,b1) = SU[i1,i1,T1,:]

            if i2 == i2old &&
                coupling_rules(a1,i3,a2) == 1 &&
                coupling_rules(b1,i3,b2) == 1

                T2 = Uf[a1,b1,i3,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,T1,T2,U3,F)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i2,i1,a1,a2,i3,m3) * sixjr(i2,i1,b1,b2,i3,n3) *
                    R_matrix(i2,i3,m3) * conj(R_matrix(i2,i3,n3))

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_3_vec(state_6p::Array{Complex{Float64},1})

    #How to combine punctures into indices?

    Amp_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

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
                    sixjr(i2,i3,m3,i1,a2,a1) * sixjr(i2,i3,n3,i1,b2,b1)

            end

        end

    end

    return chop.(Amp_new)

end




function transform_6p_state_other_1_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i2,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,m3,n3) = SU[i2,i2,U1,:]

        (i1,a2,b2) = SU[m3,n3,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i6,i5) = SF[a3,b3,F,:]

        for T1 in 1:RangeU[i1,i1]

            (i2old,a1,b1) = SU[i1,i1,T1,:]

            if  i2 == i2old &&
                coupling_rules(a1,a2,i3) == 1 &&
                coupling_rules(b1,b2,i3) == 1

                T2 = Uf[a1,b1,i3,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,T1,T2,U3,F)

                    Ampl_new[index] += state_6p_other[index_old] *
                        sixjr(i3,a2,a1,i1,i2,m3) * sixjr(i3,b2,b1,i1,i2,n3) *
                        conj(R_matrix(i2,i3,m3)) * R_matrix(i2,i3,n3)

            end

        end

    end

    return chop.(Ampl_new)

end




function transform_6p_state_other_2_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

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
                    sixjr(i3,i2,m3,a2,i1,a1) * sixjr(i3,i2,n3,b2,i1,b1)

            end

        end

    end

    return chop.(Ampl_new)

end
