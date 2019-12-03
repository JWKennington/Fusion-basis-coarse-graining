#Define the lengths for one final 6p state. Indices are the following structure:


function length_i_U_i_U_X_F()

    len = zeros(Int,k+1)

    for i in 1:k+1

        len[i] = sum(len_U_i_U_X_F[i,i,:])

    end

    return len

end

len_i_U_i_U_X_F = length_i_U_i_U_X_F()

#println(len_i_U_i_U_X_F)

#println(sum(len_i_U_i_U_X_F[:]))

#println(configs_6p)
