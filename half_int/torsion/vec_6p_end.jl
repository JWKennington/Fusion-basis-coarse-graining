#Define the lengths for one final 6p state. Indices are the following structure:
#s5,K5,U1,s6,K6,U2,X,U3,E

#println(len_U_s_K_U_X_U_E) already defined

function length_ij_U_ij_U_X_U()

    len = zeros(Int,k+1,k+1)

    for i in 1:k+1, j in 1:k+1

        len[i,j] = sum(len_U_ij_U_X_U[i,j,:])

    end

    return len

end

len_ij_U_ij_U_X_U = length_ij_U_ij_U_X_U()

#println(len_s_K_U_s_K_U_X_U_E)

#println(sum(len_ij_U_ij_U_X_U[:,:]))

#println(configs_6p)
