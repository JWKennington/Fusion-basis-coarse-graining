module Fusion

#println("Hello world")

using LinearAlgebra

include("constants.jl") #constants of quantum group
include("qgroup_def.jl") #basic functions of quantum group
include("basic_superindices.jl") #superindices for quantum group
include("vec_6p.jl") #stored length of indices for vectorizing 6p states
include("vec_10p.jl") #stored length of indices for vectorizing 10p states

include("vec_10p_SVD.jl") #stored length of indices for the svd

include("vec_8p.jl")

include("vec_8p_SVD.jl")

include("vec_6p_end.jl")


include("translate_6p.jl") #functions translating between superindices and vectorized form for 6p states.
#include("basis_trafo_6p.jl") #basis transformations 6p states
include("basis_trafo_6p_vec.jl") #basis transformation 6p states vectorized
include("basis_trafo_6p_inv_vec.jl")
include("translate_10p.jl") #functions translating between superindices and vectorized form for 10p states.
include("translate_10p_SVD.jl") #functions translating for first SVD.

include("translate_8p.jl")

include("translate_8p_SVD.jl")

include("translate_6p_end.jl")

include("basis_trafo_10p_vec.jl") #vectorized basis trafos for 10p states

include("svd_1.jl") #first SVD functions

include("basis_trafo_8p_vec.jl")

include("svd_2.jl")

include("basis_trafo_6p_end_vec.jl")


include("compute_exp.jl")



function main_vec(g::Float64)

    iterations = 15; #

    listofsv = zeros(k+1,k+1,iterations,2)

    exp_values = complex(zeros(k+1,k+1,iterations+1))

    println("Coupling constant g = ",g)

    println("Define initial state:")

    init_state = define_6p_strong_vec(g)

    exp_values[:,:,1] = chop.(compute_expectation_values(init_state))

    println("Begin iterations:")

    for iter in 1:iterations

        println()
        println()
        println("Iteration ",iter)

        init_state_other = define_6p_state_other_vec(init_state)

        init_state_1 = transform_6p_state_1_vec(init_state)

        init_state = 0.;

        init_state_2 = transform_6p_state_2_vec(init_state_1)


        init_state_1 = 0.;

        init_state_3 = transform_6p_state_3_vec(init_state_2)

        init_state_2 = 0.;

        init_state_other_1 = transform_6p_state_other_1_vec(init_state_other)

        init_state_other = 0.;

        init_state_other_2 = transform_6p_state_other_2_vec(init_state_other_1)

        init_state_other_1 = 0.;

        println()

        state_10p = glueing_6p_states_alt_vec(init_state_3, init_state_other_2)

        println()

        init_state_3 = 0.;
        init_state_other_2 = 0.;

        state_10p_1 = ten_p_trafo_1_alt_vec(state_10p)

        state_10p = 0.;

        state_10p_2 = ten_p_trafo_2_alt_vec(state_10p_1)

        state_10p_1 = 0.;

        state_10p_3 = ten_p_trafo_3_alt_vec(state_10p_2)

        state_10p_2 = 0.;

        state_10p_4 = ten_p_trafo_4_alt_vec(state_10p_3)

        state_10p_3 = 0.;

        state_10p_5 = ten_p_trafo_5_alt_vec(state_10p_4)

        state_10p_4 = 0.;

        state_10p_6 = ten_p_trafo_6_alt_vec(state_10p_5)

        state_10p_5 = 0.;

        state_10p_bSVD_alt = ten_p_trafo_bSVD_alt_vec(state_10p_6)

        state_10p_6 = 0.;

        println()
        println("Define matrix and perform SVD:")

        (embmap,newtensor,sv1) = define_matrix_for_SVD_1_vec(state_10p_bSVD_alt)

        listofsv[:,:,iter,1] = sv1

        state_8p_aSVD = define_8p_state_vec(newtensor, embmap)

        embmap = 0.;
        newtensor = 0.;

        println()

        state_8p_1 = transform_8p_state_1_vec(state_8p_aSVD)

        state_8p_aSVD = 0.;

        state_8p_2 = transform_8p_state_2_vec(state_8p_1)

        state_8p_1 = 0.;

        state_8p_3 = transform_8p_state_3_vec(state_8p_2)

        state_8p_2 = 0.;

        state_8p_4 = transform_8p_state_4_vec(state_8p_3)

        state_8p_3 = 0.;

        state_8p_5 = transform_8p_state_5_vec(state_8p_4)

        state_8p_4 = 0.;

        state_8p_6 = transform_8p_state_6_vec(state_8p_5)

        state_8p_5 = 0.;

        println()

        state_8p_bSVD = state_before_SVD_2_vec(state_8p_6)

        state_8p_6 = 0.;

        println()
        println("Define matrix and perform second SVD:")

        (embmap2,newtensor2,sv2) = define_matrix_for_SVD_2_vec(state_8p_bSVD)

        listofsv[:,:,iter,2] = sv2

        state_8p_bSVD = 0.;

        state_6p_end = define_6p_state_vec(newtensor2, embmap2)

        newtensor2 = 0.;
        embmap2 = 0.;

        println()

        state_6p_end_1 = transform_6p_state_end_1_vec(state_6p_end)

        state_6p_end = 0.;

        state_6p_end_2 = transform_6p_state_end_2_vec(state_6p_end_1)

        state_6p_end_1 = 0.;

        state_6p_end_3 = transform_6p_state_end_3_vec(state_6p_end_2)

        state_6p_end_2 = 0.;

        init_state = transform_6p_state_end_4_vec(state_6p_end_3)

        state_6p_end_3 = 0.;

        exp_values[:,:,iter+1] = chop.(compute_expectation_values(init_state))


    end

    println()
    println("Done.")
    println()

    println("Singular values:")
    println()

    for x in 1:iterations

        for i in 1:k+1, j in 1:k+1

            print(listofsv[i,j,x,1],",")

        end

        for i in 1:k+1, j in 1:k+1

            print(listofsv[i,j,x,2],",")

        end

    end

    println()
    println()
    println("Expectation values of Ribbon operators")
    println()
    println()

    for x in 1:iterations+1

        for i in 1:k+1, j in 1:k+1

            print(real(exp_values[i,j,x]),",")

        end

    end

end

@time main_vec(1.15)

end
