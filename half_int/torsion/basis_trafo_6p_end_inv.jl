function transform_6p_state_end_1_inv(state_6p::Array{Complex{Float64},7})

    Amp_new = complex(zeros(k+1,maxK,maxU,maxU,maxU,maxU,maxE))

    for s2 in 1:k+1, K2 in 1:RangeK[s2]

        (i2,j2) = S[s2,K2,:]

        for U1 in 1:RangeU[i2,j2]

            (s5,K5,M1,N1) = SU[i2,j2,U1,:]

            (i5,j5) = S[s5,K5,:]

            for U2 in 1:RangeU[M1,N1]

                (s1,K1,M5,N5) = SU[M1,N1,U2,:]

                (i1,j1) = S[s1,K1,:]

                for U3 in 1:RangeU[M5,N5]

                    (s3,K3,M3,N3) = SU[M5,N5,U3,:]

                    (i3,j3) = S[s3,K3,:]

                    for U4 in 1:RangeU[M3,N3]

                        (s6,K6,i4,j4) = SU[M3,N3,U4,:]

                        (i6,j6) = S[s6,K6,:]

                        for E in 1:RangeE[i4,j4]

                            s4 = SE[i4,j4,E]

                            for T1 in 1:RangeU[i5,j5]

                                (s1old,K1old,a1,b1) = SU[i5,j5,T1,:]

                                if s1 == s1old &&
                                    K1 == K1old &&
                                    coupling_rules(a1,i2,M5) == 1 &&
                                    coupling_rules(b1,j2,N5) == 1

                                    T2 = Uf[a1,b1,s2,K2,M5,N5]

                                    Amp_new[s2,K2,U1,U2,U3,U4,E] +=
                                    state_6p[s5,K5,T1,T2,U3,U4,E] *
                                    sixjr(i5,i2,M1,M5,i1,a1) * sixjr(j5,j2,N1,N5,j1,b1) #*
                                    #R_matrix(i4,j4,s4) * conj(R_matrix(i6,j6,s6))

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return Amp_new

end



function transform_6p_state_end_2_inv(state_6p::Array{Complex{Float64},9})

    Amp_new = complex(zeros(k+1,maxK,maxU,maxU,maxU,maxU,maxE))

    for s5 in 1:k+1, K5 in 1:RangeK[s5]

        (i5,j5) = S[s5,K5,:]

        for U1 in 1:RangeU[i5,j5]

            (s1,K1,a1,b1) = SU[i5,j5,U1,:]

            (i1,j1) = S[s1,K1,:]

            for U2 in 1:RangeU[a1,b1]

                (s2,K2,M5,N5) = SU[a1,b1,U2,:]

                (i2,j2) = S[s2,K2,:]

                for U3 in 1:RangeU[M5,N5]

                    (s3,K3,M3,N3) = SU[M5,N5,U3,:]

                    (i3,j3) = S[s3,K3,:]

                    for U4 in 1:RangeU[M3,N3]

                        (s6,K6,i4,j4) = SU[M3,N3,U4,:]

                        (i6,j6) = S[s6,K6,:]

                        for E in 1:RangeE[i4,j4]

                            s4 = SE[i4,j4,E]

                            for T3 in 1:RangeU[i6,j6]

                                (s3old,K3old,a4,b4) = SU[i6,j6,T3,:]

                                if s3 == s3old &&
                                    K3 == K3old &&
                                    coupling_rules(a4,M5,i4) == 1 &&
                                    coupling_rules(b4,N5,j4) == 1

                                    X = Xf[M5,N5,a4,b4,i4,j4]

                                    Amp_new[s5,K5,U1,U2,U3,U4,E] +=
                                    state_6p[s5,K5,U1,U2,s6,K6,T3,X,E] *
                                    sixjr(i3,M5,M3,i4,i6,a4) * sixjr(j3,N5,N3,j4,j6,b4) *
                                    R_matrix(i6,i3,a4) * conj(R_matrix(j6,j3,b4))

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return Amp_new

end



function transform_6p_state_end_3_inv(state_6p::Array{Complex{Float64},9})

    Amp_new = complex(zeros(k+1,maxK,maxU,maxU,k+1,maxK,maxU,maxX,maxE))

    for s5 in 1:k+1, K5 in 1:RangeK[s5]

        (i5,j5) = S[s5,K5,:]

        for U1 in 1:RangeU[i5,j5]

            (s1,K1,a1,b1) = SU[i5,j5,U1,:]

            (i1,j1) = S[s1,K1,:]

            for U2 in 1:RangeU[a1,b1]

                (s2,K2,M5,N5) = SU[a1,b1,U2,:]

                (i2,j2) = S[s2,K2,:]

                for s6 in 1:k+1, K6 in 1:RangeK[s6]

                    (i6,j6) = S[s6,K6,:]

                    for U3 in 1:RangeU[i6,j6]

                        (s3,K3,a4,b4) = SU[i6,j6,U3,:]

                        (i3,j3) = S[s3,K3,:]

                        for X in 1:RangeX[M5,N5,a4,b4]

                            (i4,j4) = SX[M5,N5,a4,b4,X,:]

                            for E in 1:RangeE[i4,j4]

                                s4 = SE[i4,j4,E]

                                for Xold in 1:RangeX[a1,b1,a4,b4]

                                    (a3,b3) = SX[a1,b1,a4,b4,Xold,:]

                                    if coupling_rules(i2,a3,i4) == 1 &&
                                        coupling_rules(j2,b3,j4) == 1

                                        T3 = Uf[a3,b3,s2,K2,i4,j4]

                                        Amp_new[s5,K5,U1,U2,s6,K6,U3,X,E] +=
                                        state_6p[s5,K5,U1,s6,K6,U3,Xold,T3,E] *
                                        sixjr(a1,i2,M5,i4,a4,a3) * sixjr(b1,j2,N5,j4,b4,b3)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return Amp_new

end



function transform_6p_state_end_4_inv(state_6p::Array{Complex{Float64},7})

    Amp_new = complex(zeros(k+1,maxK,maxU,k+1,maxK,maxU,maxX,maxU,maxE))

    for s5 in 1:k+1, K5 in 1:RangeK[s5]

        (i5,j5) = S[s5,K5,:]

        for U1 in 1:RangeU[i5,j5]

            (s1,K1,a1,b1) = SU[i5,j5,U1,:]

            (i1,j1) = S[s1,K1,:]

            for s6 in 1:k+1, K6 in 1:RangeK[s6]

                (i6,j6) = S[s6,K6,:]

                for U2 in 1:RangeU[i6,j6]

                    (s3,K3,a4,b4) = SU[i6,j6,U2,:]

                    (i3,j3) = S[s3,K3,:]

                    for X in 1:RangeX[a1,b1,a4,b4]

                        (a3,b3) = SX[a1,b1,a4,b4,X,:]

                        for U3 in 1:RangeU[a3,b3]

                            (s2,K2,i4,j4) = SU[a3,b3,U3,:]

                            (i2,j2) = S[s2,K2,:]

                            for E in 1:RangeE[i4,j4]

                                s4 = SE[i4,j4,E]

                                for T2 in 1:RangeU[a1,b1]

                                    (s6old,K6old,a2,b2) = SU[a1,b1,T2,:]

                                    if s6 == s6old &&
                                        K6 == K6old &&
                                        coupling_rules(a2,i3,a3) == 1 &&
                                        coupling_rules(b2,j3,b3) == 1

                                        T3 = Uf[a2,b2,s3,K3,a3,b3]

                                        Amp_new[s5,K5,U1,s6,K6,U2,X,U3,E] +=
                                        state_6p[s5,K5,U1,T2,T3,U3,E] *
                                        sixjr(a1,a3,a4,i3,i6,a2) * sixjr(b1,b3,b4,j3,j6,b2)

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return Amp_new

end
