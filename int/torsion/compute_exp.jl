function compute_expectation_values(state_6p::Array{Complex{Float64},1})

    #The idea is to compute expectation values of Ribbon operators after an iteration of the algorithm

    #To this end, we need to identify / glue opposite punctures.

    klev = div(p,2) - 2

    Z = complex(0.)
    exp = complex(zeros(klev+1,klev+1))

    for index in 1:configs_6p

        (i1,j1,U1,U2,U3,U4) = translate_6p_vec_to_ind(index)

        (i2,j2,a1,b1) = SU[i1,j1,U1,:]

        (i3,j3,a2,b2) = SU[a1,b1,U2,:]

        (i4,j4,a3,b3) = SU[a2,b2,U3,:]

        (i5,j5,i6,j6) = SU[a3,b3,U4,:]

        if i1 == i3 && j1 == j3 &&
            i2 == i4 && j2 == j4 &&
            i5 == i6 && j5 == j6

            Z += state_6p[index] * D^6 / v(i6) / v(j6) /
             v(i1) / v(j1) /
             v(i2) / v(j2)

            for i in 1:klev+1, j in 1:klev+1

                exp[i,j] += state_6p[index] * D^6 / v(i6) / v(j6) /
                v(i1) / v(j1) /
                v(i2) / v(j2) * Ribbon_exp(i,j,i1,j1)

            end

        end

    end

    return exp / Z

end
