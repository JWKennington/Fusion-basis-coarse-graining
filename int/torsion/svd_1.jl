function define_matrix_for_SVD_1_vec(state_10p::Array{Complex{Float64},1})
    mapU = complex(zeros(k+1,k+1,maxind1))

    maps = zeros(k+1,k+1)

    mapV = complex(zeros(k+1,k+1,maxind2))

    embmap = complex(zeros(k+1,k+1,maxind1))

    temp_phase = complex(0.)

    new_ten = complex(zeros(k+1,k+1,maxind2))

    for x4 in 1:k+1, y4 in 1:k+1

        mat = complex(zeros(configs_10p_ind1_svd[x4,y4],configs_10p_ind2_svd[x4,y4]))

        for index1 in 1:configs_10p_ind1_svd[x4,y4]

            T1 = translate_10p_SVD_ind1_vec_to_ind(index1,x4,y4)

            (i4,j4,k4,l4) = SU[x4,y4,T1,:]

            #Threads.@threads for index2 in 1:configs_10p_ind2_svd[x4,y4]

            for index2 in 1:configs_10p_ind2_svd[x4,y4]

                (i1,j1,U1,i2,j2,Z2,X2,U4,U5,X4,U8) = translate_10p_SVD_ind2_vec_to_ind(index2,x4,y4)

                (i3,j3,a1,b1) = SU[i1,j1,U1,:]

                (k2,l2,x2,y2) = SU[i2,j2,Z2,:]

                (a3,b3) = SX[a1,b1,x2,y2,X2,:]

                (i5,j5,M2,N2) = SU[a3,b3,U4,:]

                (k6,l6,c2,d2) = SU[M2,N2,U5,:]

                (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                (k3,l3,k1,l1) = SU[c1,d1,U8,:]

                Z4 = Uf[i4,j4,k4,l4,x4,y4]

                index_old = translate_10p_ind_to_vec_bSVD(i1,j1,U1,i2,j2,Z2,X2,U4,U5,i4,j4,Z4,X4,U8)

                mat[index1,index2] = state_10p[index_old]

            end

        end

        decomp = svd(mat)

        println("x4=",x4,", y4=",y4)

        println(decomp.S)

        mapU[x4,y4,1:configs_10p_ind1_svd[x4,y4]] = chop.(decomp.U[1:configs_10p_ind1_svd[x4,y4],1])
        maps[x4,y4] = chop.(decomp.S[1])
        mapV[x4,y4,1:configs_10p_ind2_svd[x4,y4]] = chop.(transpose(conj(decomp.V))[1,1:configs_10p_ind2_svd[x4,y4]])

        count = 0

        temp_index1 = 1

        while count < 1

            if abs(mapU[x4,y4,temp_index1]) > TOLERANCE

                count += 1

                T1 = translate_10p_SVD_ind1_vec_to_ind(temp_index1,x4,y4)

                (i4,j4,k4,l4) = SU[x4,y4,T1,:]

                temp_phase = conj(mapU[x4,y4,temp_index1] / abs(mapU[x4,y4,temp_index1])) #*
                    #phase_factor_emb_full(u4,x4,y4,s4,i4,j4,t4,k4,l4))

            end

        temp_index1 += 1

        end

        embmap[x4,y4,:] = mapU[x4,y4,:] * temp_phase


        new_ten[x4,y4,:] =
            mapV[x4,y4,:] * maps[x4,y4] / maps[1,1] *
            conj(temp_phase)


    end

    println()
    println("List of first singular values from various blocks:")

    println(maps / maps[1,1])
    #println(mapU)
    #println(mapV)


    #return(mapU,maps,mapV)
    return(chop.(embmap),chop.(new_ten),chop.(maps / maps[1,1]))

end



function define_8p_state_vec(new_ten::Array{Complex{Float64},3},emb::Array{Complex{Float64},3})

    Amp_new = complex(zeros(configs_8p))

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (i1,j1,U1,U2,U3,U4,U5,U6) = translate_8p_vec_to_ind_1(index)

        (i3,j3,a1,b1) = SU[i1,j1,U1,:]

        (x2,y2,a3,b3) = SU[a1,b1,U2,:]

        (i5,j5,M2,N2) = SU[a3,b3,U3,:]

        (k6,l6,c2,d2) = SU[M2,N2,U4,:]

        (x4,y4,c1,d1) = SU[c2,d2,U5,:]

        (k3,l3,k1,l1) = SU[c1,d1,U6,:]

        for i2 in 1:k+1, j2 in 1:k+1

            for Z2 in 1:RangeU[i2,j2]

                (k2,l2,x2old,y2old) = SU[i2,j2,Z2,:]

                if x2 == x2old &&
                    y2 == y2old

                    X2 = Xf[a1,b1,x2,y2,a3,b3]

                    T1 = Uf[x2,y2,i2,j2,k2,l2]

                    X4 = Xf[c2,d2,x4,y4,c1,d1]

                    index1_old = translate_10p_SVD_ind1_ind_to_vec(T1,x2,y2)

                    index2_old = translate_10p_SVD_ind2_ind_to_vec(i1,j1,U1,i2,j2,Z2,X2,U3,U4,X4,U6,x4,y4)

                    Amp_new[index] +=
                        new_ten[x4,y4,index2_old] *
                        emb[x2,y2,index1_old] *
                        conj(R_matrix(x2,a1,a3)) * R_matrix(y2,b1,b3) *
                        R_matrix(x4,c1,c2) * conj(R_matrix(y4,d1,d2)) *
                        v(x4) * v(y4) * v(x2) * v(y2) 

                end

            end

        end

    end

    return chop.(Amp_new)

end
