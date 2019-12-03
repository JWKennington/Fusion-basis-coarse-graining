function define_matrix_for_SVD_2_vec(state_8p::Array{Complex{Float64},1})
    #k+1,maxK,maxU,k+1,maxK,maxU,maxE,maxX,maxU,maxU,k+1,maxK,maxU,maxE,maxX,maxU,maxE
    #mat = complex(zeros((k+1)*maxK*maxU,(k+1)*maxK*maxU*(k+1)*maxK*maxU*maxE*maxX*maxU*maxU*maxX*maxU*maxE))

    mapU = complex(zeros(k+1,k+1,maxind1_8p))

    maps = zeros(k+1,k+1)

    mapV = complex(zeros(k+1,k+1,maxind2_8p))

    embmap = complex(zeros(k+1,k+1,maxind1_8p))
    #embmap[u3,Keff3,T1,E1]

    temp_phase = complex(0.)

    new_ten = complex(zeros(k+1,k+1,maxind2_8p))
    #new_ten[u3,Keff3,u2,K2,U1,s1,K1,U2,E1p,X1,X2,U4,E]

    for x3 in 1:k+1, y3 in 1:k+1

        #define matrix
        mat = complex(zeros(configs_8p_ind1_svd[x3,y3],configs_8p_ind2_svd[x3,y3]))

        for index1 in 1:configs_8p_ind1_svd[x3,y3]

            T1 = translate_8p_SVD_ind1_vec_to_ind(index1,x3,y3)

            (i3,j3,k3,l3) = SU[x3,y3,T1,:]

            #Threads.@threads for index2 in 1:configs_8p_ind2_svd[x3,y3]

            for index2 in 1:configs_8p_ind2_svd[x3,y3]

                (x2,y2,U1,i1,j1,U2,X1,X2,U4) = translate_8p_SVD_ind2_vec_to_ind(index2,x3,y3)

                (i5,j5,M1,N1) = SU[x2,y2,U1,:]

                (k1,l1,x1,y1) = SU[i1,j1,U2,:]

                (M5,N5) = SX[M1,N1,x1,y1,X1,:]

                (M3,N3) = SX[M5,N5,x3,y3,X2,:]

                (k6,l6,x4,y4) = SU[M3,N3,U4,:]

                U3 = Uf[i3,j3,k3,l3,x3,y3]

                index_old = translate_8p_ind_to_vec_bSVD(x2,y2,U1,i1,j1,U2,X1,i3,j3,U3,X2,U4)

                mat[index1,index2] = state_8p[index_old]

            end

        end

        #decomp = svds(mat,nsv=1,ritzvec=true,tol=TOLERANCE,maxiter=1000)[1]

        #mapU[u3,Keff3,:] = chop.(decomp[:U][:])
        #maps[u3,Keff3] = chop.(decomp[:S][1])
        #mapV[u3,Keff3,:] = chop.(transpose(conj(decomp[:V][:])))

        decomp = svd(mat)

        println("x3=",x3,", y3=",y3)

        #println(decomp[1])
        println(decomp.S)

        #println(decomp[3][1,:])

        mapU[x3,y3,1:configs_8p_ind1_svd[x3,y3]] = chop.(decomp.U[1:configs_8p_ind1_svd[x3,y3],1])
        maps[x3,y3] = chop.(decomp.S[1])
        mapV[x3,y3,1:configs_8p_ind2_svd[x3,y3]] = chop.(transpose(conj(decomp.V))[1,1:configs_8p_ind2_svd[x3,y3]])


        #Can one use reshape?

        #translate the maps back to previous indices
        #Threads.@threads for index1 in 1:configs_8p_ind1_svd[x3,y3]

        count = 0

        temp_index1 = 1

        while count < 1

            if abs(mapU[x3,y3,temp_index1]) > TOLERANCE

                count += 1

                T1 = translate_10p_SVD_ind1_vec_to_ind(temp_index1,x3,y3)

                (i3,j3,k3,l3) = SU[x3,y3,T1,:]

                temp_phase = conj(mapU[x3,y3,1] / abs(mapU[x3,y3,1])) #*
                    #phase_factor_emb_full(u3,x3,y3,s3,i3,j3,t3,k3,l3))

            end

        temp_index1 += 1

        end

        embmap[x3,y3,:] = mapU[x3,y3,:] * temp_phase

        new_ten[x3,y3,:] =
            mapV[x3,y3,:] * maps[x3,y3] / maps[1,1] * conj(temp_phase)

    end

    println()
    println("List of first singular values from different blocks:")

    println(maps / maps[1,1])

    #return(mapU,maps,mapV)
    return(chop.(embmap),chop.(new_ten),chop.(maps / maps[1,1]))

end


function define_6p_state_vec(new_ten::Array{Complex{Float64},3},emb::Array{Complex{Float64},3})

    Amp_new = complex(zeros(configs_6p))

    #Threads.@threads for index in 1:configs_6p
    for index in 1:configs_6p

        (x2,y2,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i5,j5,M1,N1) = SU[x2,y2,U1,:]

        (x1,y1,M5,N5) = SU[M1,N1,U2,:]

        (x3,y3,M3,N3) = SU[M5,N5,U3,:]

        (k6,l6,x4,y4) = SU[M3,N3,U4,:]

        for i1 in 1:k+1, j1 in 1:k+1

            for Z1 in 1:RangeU[i1,j1]

                (k1,l1,x1old,y1old) = SU[i1,j1,Z1,:]

                if x1 == x1old &&
                    y1 == y1old

                    X1 = Xf[M1,N1,x1,y1,M5,N5]

                    T1 = Uf[x1,y1,i1,j1,k1,l1]

                    X2 = Xf[M5,N5,x3,y3,M3,N3]

                    index1_old = translate_8p_SVD_ind1_ind_to_vec(T1,x1,y1)

                    index2_old = translate_8p_SVD_ind2_ind_to_vec(x2,y2,U1,i1,j1,Z1,X1,X2,U4,x3,y3)

                    Amp_new[index] +=
                        new_ten[x3,y3,index2_old] *
                        emb[x1,y1,index1_old] *
                        v(x1) * v(y1) * v(x3) * v(y3) #*
                        #sqrt_minus_one_pow(u3 - x3 - y3 + 1) *
                        #Asqrt^( -(x3-1)*x3 -(y3-1)*y3 + (u3-1)*u3) *
                        #sqrt_minus_one_pow(u1 -x1 - y1 + 1) *
                        #Asqrt^( (x1-1)*x1 + (y1-1)*y1 - (u1-1)*u1)
                                        #conj(R_matrix(x2,a1,a3)) * R_matrix(y2,b1,b3) *
                                        #R_matrix(x4,c1,c2) * conj(R_matrix(y4,d1,d2))

                end

            end

        end

    end

    return chop.(Amp_new)

end
