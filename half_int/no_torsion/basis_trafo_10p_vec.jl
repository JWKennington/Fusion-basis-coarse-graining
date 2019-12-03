function glueing_6p_states_alt_vec(state_6p_1::Array{Complex{Float64},1},state_6p_2::Array{Complex{Float64},1})

    newtensor = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,U6,U7,F) = translate_10p_vec_to_ind_1(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,i6,j6) = SU[M1,N1,U4,:]

        (k6,c3,d3) = SU[i6,j6,U5,:]

        (k4,c2,d2) = SU[c3,d3,U6,:]

        (k2,c1,d1) = SU[c2,d2,U7,:]

        #(k3,k1,l1) = SU[c1,d1,U8,:]

        (k3,k1) = SF[c1,d1,F,:]

        #t1 = SE[k1,l1,E]
        #L1 = Kf[t1,k1,l1]

        T1 = Uf[k1,k1,k3,c1,d1]

        T2 = Uf[c1,d1,k2,c2,d2]

        T3 = Uf[c2,d2,k4,c3,d3]

        #T4 = Uf[c3,d3,k6,i6,i6]

        #Because we have not torsion, we do not need to sum.

        #Define Fleft and Fright, such that they agree in

        if i6 == j6

            Fleft = Ff[M1,N1,i4,i6]

            Fright = Ff[c3,d3,k6,i6]

            index_left_old = translate_6p_ind_to_vec(i1,U1,U2,U3,Fleft)
            index_right_old = translate_6p_ind_to_vec(k1,T1,T2,T3,Fright)

            newtensor[index] =
                state_6p_1[index_left_old] *
                state_6p_2[index_right_old] * D^2 / v(i6)^2

        end

    end

    return chop.(newtensor)

end

function test_phase_10(state_10p::Array{Complex{Float64},1}) #Need to change orientation of punctures

    temp = 0.

    temp2 = 0.

    for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,U6,U7,F) = translate_10p_vec_to_ind_1(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (i4,i6,j6) = SU[M1,N1,U4,:]

        (k6,c3,d3) = SU[i6,j6,U5,:]

        (k4,c2,d2) = SU[c3,d3,U6,:]

        (k2,c1,d1) = SU[c2,d2,U7,:]

        (k3,k1) = SF[c1,d1,F,:]

        temp += state_10p[index] *
                minus_one_pow(k6-i1-i2-i3-i4-i5-k1-k2-k3-k4 + 8) *
                A^( (k6-1)*k6 - (i1-1)*i1 - (i2-1)*i2 - (i3-1)*i3 - (i4-1)*i4 - (i5-1)*i5 -
                (k1-1)*k1 - (k2-1)*k2 - (k3-1)*k3 - (k4-1)*k4 )

    end

    #temp2 = sum(abs.(state_10p))

    return temp

end


function ten_p_trafo_1_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    temp = 0.

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,U6,U7,F) = translate_10p_vec_to_ind_1(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (k6,c4,d4) = SU[M1,N1,U4,:]

        (i4,c3,d3) = SU[c4,d4,U5,:]

        (k4,c2,d2) = SU[c3,d3,U6,:]

        (k2,c1,d1) = SU[c2,d2,U7,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T4 in 1:RangeU[M1,N1]

            (i4old,i6,j6) = SU[M1,N1,T4,:]

            if i4 == i4old &&
                coupling_rules(i6,k6,c3) == 1 &&
                coupling_rules(j6,k6,d3) == 1

                T5 = Uf[i6,j6,k6,c3,d3]

                index_old = translate_10p_ind_to_vec_1(i1,U1,U2,U3,T4,T5,U6,U7,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i4,M1,i6,k6,c3,c4) * sixjr(i4,N1,j6,k6,d3,d4)

            end

        end

        #if a1 != b1 || a2 != b2 || M1 != N1 || c4 != d4 || c3 != d3 || c2 != d2 ||
        #    c1 != d1

        #    temp += abs(Amp_new[index])

        #end

    end

    #println(chop(temp))

    return chop.(Amp_new)

end


function ten_p_trafo_2_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,i4,Z4,X4,U7,F) = translate_10p_vec_to_ind_2(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (k6,c4,d4) = SU[M1,N1,U4,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (c2,d2) = SX[c4,d4,x4,y4,X4,:]

        (k2,c1,d1) = SU[c2,d2,U7,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T5 in 1:RangeU[c4,d4]

            (i4old,c3,d3) = SU[c4,d4,T5,:]

            if i4 == i4old &&
                coupling_rules(c3,k4,c2) == 1 &&
                coupling_rules(d3,k4,d2) == 1

                T6 = Uf[c3,d3,k4,c2,d2]

                index_old = translate_10p_ind_to_vec_1(i1,U1,U2,U3,U4,T5,T6,U7,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i4,c4,c3,c2,k4,x4) * sixjr(i4,d4,d3,d2,k4,y4)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_3_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,i4,Z4,U7,X4,F) = translate_10p_vec_to_ind_3(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (k6,c4,d4) = SU[M1,N1,U4,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (k2,m2,n2) = SU[x4,y4,U7,:]

        (c1,d1) = SX[c4,d4,m2,n2,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        for X4old in 1:RangeX[c4,d4,x4,y4]

            (c2,d2) = SX[c4,d4,x4,y4,X4old,:]

            if coupling_rules(c2,k2,c1) == 1 &&
                coupling_rules(d2,k2,d1) == 1

                T7 = Uf[c2,d2,k2,c1,d1]

                index_old = translate_10p_ind_to_vec_2(i1,U1,U2,U3,U4,i4,Z4,X4old,T7,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(x4,c4,c2,c1,k2,m2) * sixjr(y4,d4,d2,d1,k2,n2) *
                    conj(R_matrix(k2,x4,m2)) * R_matrix(k2,y4,n2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_4_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F) = translate_10p_vec_to_ind_4(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (k6,c4,d4) = SU[M1,N1,U4,:]

        (k2,c2,d2) = SU[c4,d4,U5,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T7 in 1:RangeU[x4,y4]

            (k2old,m2,n2) = SU[x4,y4,T7,:]

            if k2 == k2old &&
                coupling_rules(c4,m2,c1) == 1 &&
                coupling_rules(d4,n2,d1) == 1

                X4old = Xf[c4,d4,m2,n2,c1,d1]

                index_old = translate_10p_ind_to_vec_3(i1,U1,U2,U3,U4,i4,Z4,T7,X4old,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(x4,k2,m2,c4,c1,c2) * sixjr(y4,k2,n2,d4,d1,d2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_5_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F) = translate_10p_vec_to_ind_4(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (i5,M1,N1) = SU[a2,b2,U3,:]

        (k2,M2,N2) = SU[M1,N1,U4,:]

        (k6,c2,d2) = SU[M2,N2,U5,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T4 in 1:RangeU[M1,N1]

            (k6old,c4,d4) = SU[M1,N1,T4,:]

            if k6 == k6old &&
                coupling_rules(c4,k2,c2) == 1 &&
                coupling_rules(d4,k2,d2) == 1

                T5 = Uf[c4,d4,k2,c2,d2]

                index_old = translate_10p_ind_to_vec_4(i1,U1,U2,U3,T4,T5,i4,Z4,X4,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(M1,k6,c4,c2,k2,M2) * sixjr(N1,k6,d4,d2,k2,N2)

            end

        end

    end

    return chop.(Amp_new)

end


function ten_p_trafo_6_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F) = translate_10p_vec_to_ind_4(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (k2,a3,b3) = SU[a2,b2,U3,:]

        (i5,M2,N2) = SU[a3,b3,U4,:]

        (k6,c2,d2) = SU[M2,N2,U5,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T3 in 1:RangeU[a2,b2]

            (i5old,M1,N1) = SU[a2,b2,T3,:]

            if i5 == i5old &&
                coupling_rules(M1,M2,k2) == 1 &&
                coupling_rules(N1,N2,k2) == 1

                T4 = Uf[M1,N1,k2,M2,N2]

                index_old = translate_10p_ind_to_vec_4(i1,U1,U2,T3,T4,U5,i4,Z4,X4,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(a2,i5,M1,M2,k2,a3) * sixjr(b2,i5,N1,N2,k2,b3)

            end

        end

    end

    return chop.(Amp_new)

end

function compare_10p_to_AL(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p
    for index in 1:configs_10p

        (i1,U1,U2,U3,U4,U5,i4,Z4,X4,F) = translate_10p_vec_to_ind_4(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (i2,a2,b2) = SU[a1,b1,U2,:]

        (k2,a3,b3) = SU[a2,b2,U3,:]

        (i5,M2,N2) = SU[a3,b3,U4,:]

        (k6,c2,d2) = SU[M2,N2,U5,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        Amp_new[index] += v(i1) * v(i2) * v(i3) * v(i4) * v(i5) *
            v(k2) * v(k3) * v(k1) * v(k4) * v(k6) / D^10 *
            minus_one_pow(k6-i1-i2-i3-i4-i5-k1-k2-k3-k4 + 8) *
            conj(A^( (k6-1)*k6 - (i1-1)*i1 - (i2-1)*i2 - (i3-1)*i3 - (i4-1)*i4 - (i5-1)*i5 -
                (k1-1)*k1 - (k2-1)*k2 - (k3-1)*k3 - (k4-1)*k4 ))

        println(Amp_new[index]," , ",state_10p[index])

    end

end


function ten_p_trafo_bSVD_alt_vec(state_10p::Array{Complex{Float64},1})

    Amp_new = complex(zeros(configs_10p))

    #Threads.@threads for index in 1:configs_10p_bSVD
    for index in 1:configs_10p
    #Threads.@threads for index in 1:configs_10p_bSVD

        (i1,U1,i2,Z2,X2,U4,U5,i4,Z4,X4,F) = translate_10p_vec_to_ind_bSVD(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (k2,x2,y2) = SU[i2,i2,Z2,:]

        #u2 = SE[x2,y2,E2]

        (a3,b3) = SX[a1,b1,x2,y2,X2,:]

        (i5,M2,N2) = SU[a3,b3,U4,:]

        (k6,c2,d2) = SU[M2,N2,U5,:]

        (k4,x4,y4) = SU[i4,i4,Z4,:]

        #u4 = SE[x4,y4,E4]

        (c1,d1) = SX[c2,d2,x4,y4,X4,:]

        (k3,k1) = SF[c1,d1,F,:]

        for T2 in 1:RangeU[a1,b1]

            (i2old,a2,b2) = SU[a1,b1,T2,:]

            if i2 == i2old &&
                coupling_rules(a2,k2,a3) == 1 &&
                coupling_rules(b2,k2,b3) == 1

                T3 = Uf[a2,b2,k2,a3,b3]

                index_old = translate_10p_ind_to_vec_4(i1,U1,T2,T3,U4,U5,i4,Z4,X4,F)

                Amp_new[index] +=
                    state_10p[index_old] *
                    sixjr(i2,a1,a2,a3,k2,x2) * sixjr(i2,b1,b2,b3,k2,y2) /
                    #v(i2) / v(k2) / v(i4) / v(k4)
                    #v(u2)^2 * v(u4)^2 / v(x2)^2 / v(y2)^2 / v(x4)^2 / v(y4)^2 / D^2
                    #v(u4) /
                    v(x4) / v(y4) / D /
                    #v(u2) /
                    v(x2) / v(y2) / D #*
                    #sqrt(v(x4)) * sqrt(v(y4)) * sqrt(v(x2)) * sqrt(v(y2))

            end

        end

    end

    return chop.(Amp_new)

end
