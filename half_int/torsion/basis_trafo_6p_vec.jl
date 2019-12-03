function define_6p_strong_vec(g::Float64)

    #g = 0.5

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i5,j5,i6,j6) = SU[a3,b3,U4,:]

        if i1 == j1 && i2 == j2 && i3 == j3 && i4 == j4 && i5 == j5 && i6 == j6 &&
            a1 == b1 && a2 == b2 && a3 == b3

            Amp_new[index] = v(i1) * v(i2) * v(i3) * v(i4) * v(i5) * v(i6) / D^6 *
                heat_ker(g,i1) * heat_ker(g,i2) * heat_ker(g,i3) * heat_ker(g,i4) * heat_ker(g,i5) * heat_ker(g,i6)

        end

    end

    return chop.(Amp_new/Amp_new[1])

end


#While we are at it, change the indices such that 'E' describes puncture 5 instead of 6.
#That way it is more easy to identify / sum over them in the glueing step of the algorithm.



function define_6p_state_other_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i6,j6,i5,j5) = SU[a3,b3,U4,:]

        T4 = Uf[a3,b3,i5,j5,i6,j6]

        index_old = translate_6p_ind_to_vec(i1,j1,U1,U2,U3,T4)

        Amp_new[index] = state_6p[index_old] *
            conj(R_matrix(i2,i1,a1)) * R_matrix(j2,j1,b1) *
            conj(R_matrix(i3,a1,a2)) * R_matrix(j3,b1,b2) *
            conj(R_matrix(i4,a2,a3)) * R_matrix(j4,b2,b3)

    end

    return chop.(Amp_new)

end



function transform_6p_state_1_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        for T3 in 1:RangeU[a2,b2]

            (i4old,j4old,a3,b3) = SU[a2,b2,T3,:]

            if  i4 == i4old &&
                j4 == j4old &&
                coupling_rules(a3,i5,i6) == 1 &&
                coupling_rules(b3,j5,j6) == 1

                T4 = Uf[a3,b3,i5,j5,i6,j6]

                index_old = translate_6p_ind_to_vec(i1,j1,U1,U2,T3,T4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i4,a2,a3,i5,i6,M1) * sixjr(j4,b2,b3,j5,j6,N1)

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_1_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i5,j5,i6,j6) = SU[a3,b3,U4,:]

        for T3 in 1:RangeU[a2,b2]

            (i5old,j5old,M1,N1) = SU[a2,b2,T3,:]

            if i5 == i5old &&
                j5 == j5old &&
                coupling_rules(M1,i4,i6) == 1 &&
                coupling_rules(N1,j4,j6) == 1

                T4 = Uf[M1,N1,i4,j4,i6,j6]

                index_old = translate_6p_ind_to_vec(i1,j1,U1,U2,T3,T4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i4,a2,a3,i5,i6,M1) * sixjr(j4,b2,b3,j5,j6,N1)

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_2_vec(state_6p::Array{Complex{Float64},1})

    #How to combine punctures into indices?

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i2,j2,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,m3,n3) = SU[i2,j2,U1,:]

        (i1,j1,a2,b2) = SU[m3,n3,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        for T1 in 1:RangeU[i1,j1]

            (i2old,j2old,a1,b1) = SU[i1,j1,T1,:]

            if i2 == i2old &&
                j2 == j2old &&
                coupling_rules(a1,i3,a2) == 1 &&
                coupling_rules(b1,j3,b2) == 1

                T2 = Uf[a1,b1,i3,j3,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,j1,T1,T2,U3,U4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i2,i1,a1,a2,i3,m3) * sixjr(j2,j1,b1,b2,j3,n3) *
                    R_matrix(i2,i3,m3) * conj(R_matrix(j2,j3,n3))

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_2_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        for T1 in 1:RangeU[i2,j2]

            (i3old,j3old,m3,n3) = SU[i2,j2,T1,:]

            if i3 == i3old &&
                j3 == j3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(j1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,j1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,j2,T1,T2,U3,U4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i2,i1,a1,a2,i3,m3) * sixjr(j2,j1,b1,b2,j3,n3) *
                    conj(R_matrix(i2,i3,m3)) * R_matrix(j2,j3,n3)

            end

        end

    end

    return chop.(Amp_new)

end


function transform_6p_state_3_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        for T1 in 1:RangeU[i2,j2]

            (i3old,j3old,m3,n3) = SU[i2,j2,T1,:]

            if i3 == i3old &&
                j3 == j3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(j1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,j1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,j2,T1,T2,U3,U4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i2,i3,m3,i1,a2,a1) * sixjr(j2,j3,n3,j1,b2,b1)

            end

        end

    end

    return chop.(Amp_new)

end


function inverse_tranform_6p_3_inv_vec(state_6p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i2,j2,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,m3,n3) = SU[i2,j2,U1,:]

        (i1,j1,a2,b2) = SU[m3,n3,U2,:]

        (i5,j5,M1,N1) = SU[a2,b2,U3,:]

        (i4,j4,i6,j6) = SU[M1,N1,U4,:]

        for T1 in 1:RangeU[i1,j1]

            (i3old,j3old,a1,b1) = SU[i1,j1,T1,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(a1,i2,a2) == 1 &&
                coupling_rules(b1,j2,b2) == 1

                T2 = Uf[a1,b1,i2,j2,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,j1,T1,T2,U3,U4)

                Amp_new[index] += state_6p[index_old] *
                    sixjr(i3,i1,a1,a2,i2,m3) * sixjr(j3,j1,b1,b2,j2,n3)

            end

        end

    end

    return chop.(Amp_new)

end



function transform_6p_state_other_1_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i2,j2,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,m3,n3) = SU[i2,j2,U1,:]

        (i1,j1,a2,b2) = SU[m3,n3,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i6,j6,i5,j5) = SU[a3,b3,U4,:]

        for T1 in 1:RangeU[i1,j1]

            (i2old,j2old,a1,b1) = SU[i1,j1,T1,:]

            if  i2 == i2old &&
                j2 == j2old &&
                coupling_rules(a1,a2,i3) == 1 &&
                coupling_rules(b1,b2,j3) == 1

                T2 = Uf[a1,b1,i3,j3,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,j1,T1,T2,U3,U4)

                    Ampl_new[index] += state_6p_other[index_old] *
                        sixjr(i3,a2,a1,i1,i2,m3) * sixjr(j3,b2,b1,j1,j2,n3) *
                        conj(R_matrix(i2,i3,m3)) * R_matrix(j2,j3,n3)

            end

        end

    end

    return chop.(Ampl_new)

end



function transform_6p_state_other_1_inv_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i6,j6,i5,j5) = SU[a3,b3,U4,:]

        for T1 in 1:RangeU[i2,j2]

            (i3old,j3old,m3,n3) = SU[i2,j2,T1,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(j1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,j1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,j2,T1,T2,U3,U4)

                Ampl_new[index] += state_6p_other[index_old] *
                    sixjr(i3,a2,a1,i1,i2,m3) * sixjr(j3,b2,b1,j1,j2,n3) *
                    R_matrix(i2,i3,m3) * conj(R_matrix(j2,j3,n3))

            end

        end

    end

    return chop.(Ampl_new)

end



function transform_6p_state_other_2_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (i2,j2,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i6,j6,i5,j5) = SU[a3,b3,U4,:]

        for T1 in 1:RangeU[i2,j2]

            (i3old,j3old,m3,n3) = SU[i2,j2,T1,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(i1,m3,a2) == 1 &&
                coupling_rules(j1,n3,b2) == 1

                T2 = Uf[m3,n3,i1,j1,a2,b2]

                index_old = translate_6p_ind_to_vec(i2,j2,T1,T2,U3,U4)

                Ampl_new[index] += state_6p_other[index_old] *
                    sixjr(i3,i2,m3,a2,i1,a1) * sixjr(j3,j2,n3,b2,j1,b1)

            end

        end

    end

    return chop.(Ampl_new)

end



function transform_6p_state_other_2_inv_vec(state_6p_other::Array{Complex{Float64},1})

    Ampl_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (i2,j2,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i3,j3,m3,n3) = SU[i2,j2,U1,:]

        (i1,j1,a2,b2) = SU[m3,n3,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i6,j6,i5,j5) = SU[a3,b3,U4,:]

        for T1 in 1:RangeU[i1,j1]

            (i3old,j3old,a1,b1) = SU[i1,j1,T1,:]

            if  i3 == i3old &&
                j3 == j3old &&
                coupling_rules(i2,a1,a2) == 1 &&
                coupling_rules(j2,b1,b2) == 1

                T2 = Uf[a1,b1,i2,j2,a2,b2]

                index_old = translate_6p_ind_to_vec(i1,j1,T1,T2,U3,U4)

                Ampl_new[index] += state_6p_other[index_old] *
                    sixjr(i3,i2,m3,a2,i1,a1) * sixjr(j3,j2,n3,b2,j1,b1)

            end

        end

    end

    return chop.(Ampl_new)

end
