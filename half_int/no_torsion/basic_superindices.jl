function compute_RangeK_S_Kf()

    RangeK = zeros(Int64,k+1)
    S = zeros(Int64,k+1,(k+1)^2,2)
    Kf = zeros(Int64,k+1,k+1,k+1)

    # RangeK[:] = typemax(Int)
    # S[:,:,:] = typemax(Int)
    # Kf[:,:,:] = typemax(Int)
    for j in 1:k+1
        K = 1
        for j1 in 1:k+1, j2 in 1:k+1
            Kf[j,j1,j2] = 1
            if coupling_rules(j1,j2,j) == 1

                RangeK[j] = K

                S[j,K,1] = j1
                S[j,K,2] = j2

                Kf[j,j1,j2] = K

                K += 1
            end
        end
        #@assert K <= z
    end

    #for j in 0:k
        #println("RangeK[$j]=$(RangeK[j])")
    #end

    return (RangeK,S,Kf)

end


const (RangeK,S,Kf)=compute_RangeK_S_Kf();

function max_Range_K(RangeK::Array{Int64,1})

    max = 0

    for i in 1:length(RangeK)

        if RangeK[i] > max

            max = RangeK[i]

        end

    end

    return max

end


const maxK = max_Range_K(RangeK);


function compute_RangeB_Jf()
    # RangeB[:,:,:] = typemax(Int)
    # Jf[:,:,:,:] = typemax(Int)

    RangeB = zeros(Int64,k+1,maxK,maxK)
    Jf = zeros(Int64,k+1,maxK,maxK,(k+1))
    Bf = zeros(Int64,k+1,maxK,maxK,k+1)

    for j5 in 1:k+1
        for K in 1:RangeK[j5], L in 1:RangeK[j5]
            B = 1
            for j6 in 1:k+1
                if (
                        Kf[j6, S[j5,K,2], S[j5,L,1]] >= 1 &&
                        Kf[j6, S[j5,L,2], S[j5,K,1]] >= 1)
                    Jf[j5,K,L,B] = j6
                    RangeB[j5,K,L] = B
                    Bf[j5,K,L,j6] = B
                    B += 1
                end
            end
            #@assert B <= k + 1
            #if B > 0
                #RangeB[j5,K,L] = B
                ##println("RangeB[$j5,$K,$L]=$(RangeB[j5,K,L])")
            #end
        end
    end

    return (RangeB,Jf,Bf)

end


const (RangeB,Jf,Bf) = compute_RangeB_Jf();


function max_Range_B()

    max = 0

    for j5 in 1:k+1

        for K in 1:RangeK[j5], L in 1:RangeK[j5]

            if RangeB[j5,K,L] > max

                max = RangeB[j5,K,L]

            end

        end

    end

    return max

end


const maxB = max_Range_B();


function compute_Fg()

    Fg = zeros(k+1,maxK,maxK,maxB)

    for j in 1:k+1
        for K in 1:RangeK[j], L in 1:RangeK[j]
            for b in 1:RangeB[j,K,L]
                Fg[j,K,L,b] =
                    chop(sixjr(
                        S[j,L,2], S[j,K,1], Jf[j,K,L,b], S[j,K,2], S[j,L,1], j))
                #println("Fg[$j,$K,$L,$b]=$(Fg[j,K,L,b])")
            end
        end
    end

    return Fg

end


const Fg = compute_Fg();


function R_matrix(i::Int,j::Int,k::Int)

    chop(minus_one_pow(k - i - j + 1) * A^(k^2 - k - i^2 + i - j^2 + j))

end


# compute maximal length of this index first

function compute_univ_index_max()

    maxU = 0

    for a1 in 1:k+1, b1 in 1:k+1

        U = 1

        for i in 1:k+1, a2 in 1:k+1, b2 in 1:k+1

            if coupling_rules(a1,i,a2) == 1 &&
                coupling_rules(b1,i,b2) == 1

                if U > maxU

                    maxU = U

                end

                U += 1

            end

        end

    end

    return maxU

end


const maxU = compute_univ_index_max();


function compute_univ_index()

    RangeU = zeros(Int,k+1,k+1)
    SU = zeros(Int,k+1,k+1,maxU,3)
    Uf = zeros(Int,k+1,k+1,k+1,k+1,k+1)

    for a1 in 1:k+1, b1 in 1:k+1

        U = 1

        for i in 1:k+1, a2 in 1:k+1, b2 in 1:k+1

            if coupling_rules(a1,i,a2) == 1 &&
                coupling_rules(b1,i,b2) == 1

                RangeU[a1,b1] = U

                SU[a1,b1,U,1] = i
                SU[a1,b1,U,2] = a2
                SU[a1,b1,U,3] = b2

                Uf[a1,b1,i,a2,b2] = U

                U += 1

            end

        end

    end

    return (RangeU,SU,Uf)

end


const (RangeU,SU,Uf) = compute_univ_index();


function compute_max_length_close()

    maxE = 0

    for a in 1:k+1, b in 1:k+1

        E = 1

        for s in 1:k+1

            if coupling_rules(a,b,s) == 1

                if E > maxE

                    maxE = E

                end

                E += 1

            end

        end

    end

    return maxE

end


const maxE = compute_max_length_close()


function compute_ind_close()

    RangeE = zeros(Int,k+1,k+1)
    SE = zeros(Int,k+1,k+1,maxE)
    Ef = zeros(Int,k+1,k+1,k+1)

    for a in 1:k+1, b in 1:k+1

        E = 1

        for s in 1:k+1

            if coupling_rules(a,b,s) == 1

                RangeE[a,b] = E

                SE[a,b,E] = s

                Ef[a,b,s] = E

                E += 1

            end

        end

    end

    return (RangeE,SE,Ef)

end


const (RangeE,SE,Ef) = compute_ind_close();


function compute_max_three_valent_length()

    maxX = 0

    for a1 in 1:k+1, b1 in 1:k+1, m3 in 1:k+1, n3 in 1:k+1

        X = 1

        for a2 in 1:k+1, b2 in 1:k+1

            if  coupling_rules(a1,m3,a2) == 1
                coupling_rules(b1,n3,b2) == 1

                if X > maxX

                    maxX = X

                end

                X += 1

            end

        end

    end

    return maxX

end


maxX = compute_max_three_valent_length();


function compute_max_three_index()

    RangeX = zeros(Int,k+1,k+1,k+1,k+1)
    SX = zeros(Int,k+1,k+1,k+1,k+1,maxX,2)
    Xf = zeros(Int,k+1,k+1,k+1,k+1,k+1,k+1)

    for a1 in 1:k+1, b1 in 1:k+1, m3 in 1:k+1, n3 in 1:k+1

        X = 1

        for a2 in 1:k+1, b2 in 1:k+1

            if  coupling_rules(a1,m3,a2) == 1 &&
                coupling_rules(b1,n3,b2) == 1

                RangeX[a1,b1,m3,n3] = X

                SX[a1,b1,m3,n3,X,1] = a2
                SX[a1,b1,m3,n3,X,2] = b2

                Xf[a1,b1,m3,n3,a2,b2] = X

                X += 1

            end

        end

    end

    return (RangeX,SX,Xf)

end


const (RangeX,SX,Xf) = compute_max_three_index();


#Compute new index for closing two punctures.

function compute_univ_close_index_max()

    maxF = 0

    for a1 in 1:k+1, b1 in 1:k+1

        F = 1

        for i in 1:k+1, i2 in 1:k+1

            if coupling_rules(a1,i,i2) == 1 &&
                coupling_rules(b1,i,i2) == 1

                if F > maxF

                    maxF = F

                end

                F += 1

            end

        end

    end

    return maxU

end


const maxF = compute_univ_close_index_max();


function compute_univ_close_index()

    RangeF = zeros(Int,k+1,k+1)
    SF = zeros(Int,k+1,k+1,maxU,2)
    Ff = zeros(Int,k+1,k+1,k+1,k+1)

    for a1 in 1:k+1, b1 in 1:k+1

        F = 1

        for i in 1:k+1, i2 in 1:k+1

            if coupling_rules(a1,i,i2) == 1 &&
                coupling_rules(b1,i,i2) == 1

                RangeF[a1,b1] = F

                SF[a1,b1,F,1] = i
                SF[a1,b1,F,2] = i2

                Ff[a1,b1,i,i2] = F

                F += 1

            end

        end

    end

    return (RangeF,SF,Ff)

end


const (RangeF,SF,Ff) = compute_univ_close_index();
