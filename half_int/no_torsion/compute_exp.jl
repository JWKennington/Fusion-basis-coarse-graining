function compute_expectation_values(state_6p::Array{Complex{Float64},1})

    #The idea is to compute expectation values of Ribbon operators after an iteration of the algorithm

    #To this end, we need to identify / glue opposite punctures.

    Z = complex(0.)
    exp = complex(zeros(k+1,k+1))

    for index in 1:configs_6p
    #Threads.@threads for index in 1:configs_6p

        (i1,U1,U2,U3,F) = translate_6p_vec_to_ind(index)

        (i2,a1,b1) = SU[i1,i1,U1,:]

        (i3,a2,b2) = SU[a1,b1,U2,:]

        (i4,a3,b3) = SU[a2,b2,U3,:]

        (i5,i6) = SF[a3,b3,F,:]

        if i1 == i3 && i2 == i4 && i5 == i6

            Z += state_6p[index] * D^6 / v(i6)^2 / v(i1)^2 / v(i2)^2

            for i in 1:k+1, j in 1:k+1

                exp[i,j] += state_6p[index] * D^6 / v(i6)^2 / v(i1)^2 / v(i2)^2 *
                    Ribbon_exp(i,j,i1,i1)

            end

        end

    end

    return exp / Z

end
