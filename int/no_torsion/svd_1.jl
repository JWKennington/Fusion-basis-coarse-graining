function define_matrix_for_SVD_1_vec(state_10p::Array{Complex{Float64},1})

    mapU = complex(zeros(k+1,k+1,maxind1))

    maps = zeros(k+1,k+1)

    mapV = complex(zeros(k+1,k+1,maxind2))

    embmap = complex(zeros(k+1,k+1,maxind1))

    temp_phase = complex(0.)

    new_ten = complex(zeros(k+1,k+1,maxind2))

    for x4 in 1:k+1, y4 in 1:k+1

        if configs_10p_ind1_svd[x4,y4] + configs_10p_ind2_svd[x4,y4] > 0

            #define matrix
            mat = complex(zeros(configs_10p_ind1_svd[x4,y4],configs_10p_ind2_svd[x4,y4]))

            for index1 in 1:configs_10p_ind1_svd[x4,y4]

                #(T1,E1) = translate_10p_SVD_ind1_vec_to_ind(index1,x4,y4)

                F1 = index1

                (i4,k4) = SF[x4,y4,F1,:]

                #Threads.@threads for index2 in 1:configs_10p_ind2_svd[x4,y4]

                for index2 in 1:configs_10p_ind2_svd[x4,y4]
                    #Threads.@threads for index2 in 1:configs_10p_ind2_svd[x4,y4]

                    (i1,U1,i2,Z2,X2,U4,U5,X4,F) = translate_10p_SVD_ind2_vec_to_ind(index2,x4,y4)

                    (i3,a1,b1) = SU[i1,i1,U1,:]

                    (k2,x2,y2) = SU[i2,i2,Z2,:]

                    #u2 = SE[x2,y2,E2p]

                    (a3,b3) = SX[a1,b1,x2,y2,X2,:]

                    (i5,M2,N2) = SU[a3,b3,U4,:]

                    (k6,c2,d2) = SU[M2,N2,U5,:]

                    (c1,d1) = SX[c2,d2,x4,y4,X4,:]

                    (k3,k1) = SF[c1,d1,F,:]

                    Z4 = Uf[i4,i4,k4,x4,y4]

                    #E4p = Ef[x4,y4,u4]

                    index_old = translate_10p_ind_to_vec_bSVD(i1,U1,i2,Z2,X2,U4,U5,i4,Z4,X4,F)

                    mat[index1,index2] = state_10p[index_old]

                end

            end

            decomp = svd(mat)

            println("x4=",x4,", y4=",y4)

            #println(decomp[1])

            println(decomp.S)
            #println(decomp[3][1,:])

            mapU[x4,y4,1:configs_10p_ind1_svd[x4,y4]] = chop.(decomp.U[1:configs_10p_ind1_svd[x4,y4],1])
            maps[x4,y4] = chop.(decomp.S[1])
            mapV[x4,y4,1:configs_10p_ind2_svd[x4,y4]] = chop.(transpose(conj(decomp.V))[1,1:configs_10p_ind2_svd[x4,y4]])

            count = 0

            temp_index1 = 1

            while count < 1

                if abs(mapU[x4,y4,temp_index1]) > TOLERANCE

                    count += 1

                    (i4,k4) = SF[x4,y4,temp_index1,:]

                    temp_phase = conj(mapU[x4,y4,temp_index1] /
                    abs(mapU[x4,y4,temp_index1])) #*
                    #phase_factor_emb_full(u4,x4,y4,1,i4,i4,1,k4,k4)

                    #println(temp_phase,",",temp_index1)

                end

                temp_index1 += 1

            end

            embmap[x4,y4,:] = mapU[x4,y4,:] * temp_phase

            new_ten[x4,y4,:] =
            mapV[x4,y4,:] * maps[x4,y4] / maps[1,1] * conj(temp_phase)


        end

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

    temp = 0.
    temp2 = 0.

    #Threads.@threads for index in 1:configs_8p
    for index in 1:configs_8p

        (i1,U1,U2,U3,U4,U5,F) = translate_8p_vec_to_ind_1(index)

        (i3,a1,b1) = SU[i1,i1,U1,:]

        (x2,a3,b3) = SU[a1,b1,U2,:]

        (i5,M2,N2) = SU[a3,b3,U3,:]

        (k6,c2,d2) = SU[M2,N2,U4,:]

        (x4,c1,d1) = SU[c2,d2,U5,:]

        (k3,k1) = SF[c1,d1,F,:]

        for Fold in 1:RangeF[x2,x2]

            (i2,k2) = SF[x2,x2,Fold,:]

            Z2 = Uf[i2,i2,k2,x2,x2]

            #E2p = Ef[x2,x2,1]

            X2 = Xf[a1,b1,x2,x2,a3,b3]

            X4 = Xf[c2,d2,x4,x4,c1,d1]

            index1_old = Fold

            index2_old = translate_10p_SVD_ind2_ind_to_vec(i1,U1,i2,Z2,X2,U3,U4,X4,F,x4,x4)


            Amp_new[index] +=
                new_ten[x4,x4,index2_old] *
                emb[x2,x2,index1_old] *
                conj(R_matrix(x2,a1,a3)) * R_matrix(x2,b1,b3) *
                R_matrix(x4,c1,c2) * conj(R_matrix(x4,d1,d2)) *
                v(x4)^2 * v(x2)^2

        end


    end

    return chop.(Amp_new)

end
