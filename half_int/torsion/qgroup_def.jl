kronecker_delta(num1::Int, num2::Int) = Int(num1 == num2)


minus_one_pow(n::Int) = ifelse(iseven(n), 1, -1)


@inline sqrt_minus_one_pow(n::Int) = ifelse(iseven(n),
                                    ifelse(iseven(div(n,2)),1,-1),
                                    ifelse(iseven(div(n-1,2)),0+1im,0-1im ))


#function chop{T<:Real}(num::T)
#    ifelse(abs(num) >= TOLERANCE, num, T(0))
#end


function chop(num::T) where T<:Real
    ifelse(abs(num) >= TOLERANCE, num, T(0))
end


#function chop{T<:Complex}(num::T)
#    T(chop(real(num)), chop(imag(num)))
#end


function chop(num::T) where T<:Complex
    T(chop(real(num)), chop(imag(num)))
end


function q0(n::Int)
    # note this could be simplified!
    An = A^n
    chop(real((An - 1/An) / (A - 1/A)))
end


const q_values = Float64[q0(n) for n in -num_q_values:num_q_values];


q(n::Int) = @inbounds return q_values[n+num_q_values+1]


qf0(n::Int) = chop(prod([q0(i) for i in 1:n]))


const qf_values = Float64[qf0(n) for n in -num_q_values:num_q_values];


qf(n::Int) = @inbounds return qf_values[n+num_q_values+1]


w0(n::Int) = chop(sqrt(q0(n)+0im))



const w_values = Complex{Float64}[w0(n) for n in -num_q_values:num_q_values];


@inline w(n::Int) = @inbounds return w_values[n+num_q_values+1];


v0(n::Int) =  chop(sqrt_minus_one_pow(n-1) * w(n));


const v_values = Complex{Float64}[v0(n) for n in -num_q_values:num_q_values];

@inline v(n::Int) = @inbounds return v_values[n+num_q_values+1];


const D = sqrt( (k+2)/2 ) / sin(pi/(k+2))


function tria(a::Int, b::Int, c::Int)
    coupling_rules(a,b,c) * #this ensures that we do not compute factorial of half integers
    sqrt(qf(div(a + b - c - 1,2)) * qf(div(a - b + c - 1,2)) * qf(div(-a + b + c - 1,2)) /
    qf(div(a + b + c - 1,2)))
end

function coupling_rules(j1::Int,j2::Int,J::Int)

    res = 0

    if  j1 + j2 - 1 >= J &&
        j1 + J  - 1 >= j2 &&
        j2 + J  - 1 >= j1 &&
        j1 + j2 + J <= p - 1 &&
        isodd(j1 + j2 + J)

        res = 1

    end

    return res

end



function sixjr(a::Int, b::Int, e::Int, d::Int, c::Int, f::Int)

    s = 0.0

    cond1 = div(a + b + e -3,2)# a + b + e
    cond2 = div(a + c + f -3,2)# a + c + f
    cond3 = div(b + d + f -3,2)# b + d + f
    cond4 = div(d + c + e -3,2)# d + c + e

    cond5 = div(a + b + c + d -4,2) # a + b + c + d
    cond6 = div(a + d + e + f -4,2) # a + d + e + f
    cond7 = div(b + c + e + f -4,2) # b + c + e + f

    lowerBound = max(cond1, cond2, cond3, cond4)
    upperBound = min(cond5, cond6, cond7)

    for n in lowerBound:upperBound
        s += minus_one_pow(n) * qf(n + 1) / qf(n - cond1) /
        qf(n - cond2) / qf(n - cond3) / qf(n - cond4) /
        qf(cond5 - n) / qf(cond6 - n) / qf(cond7 - n)
    end

    s *= sqrt(q(e) * q(f)) *
    minus_one_pow(cond5) *
    tria(a, b, e) * tria(c, d, e) *
    tria(a, c, f) * tria(b, d, f)

    return s

end


#function s_matrix(i::Int,j::Int)

#    res = minus_one_pow(2*(i+j-2)) * q((2*i-1)*(2*j-1))

#    res

#end

#function Ribbon_exp(i::Int,j::Int,a::Int,b::Int)

#    res = s_matrix(i,a) * s_matrix(j,b) / v(a)^2 / v(b)^2

#    res

#end


function Wilson(g::Float64,i::Int)

    #Define the action by comparing eigenvalues of the Wilson loop operator
    #for label l=1/2

    prefac = pi / (k + 2)

    eigval = sin(prefac) * sin(2 * prefac * i) /
        sin(2 * prefac) / sin(prefac * i)

    res = exp(1 / g^2 * (eigval - 1))

    res

end

#println(div(p,2) - 2,k)

#println(Wilson(2.,1))
#println(Wilson(2.,2))
#println(Wilson(2.,3))

#Code the heat kernel action

function flatness_op(i::Int)

    res = 0.

    prefac = pi / (k + 2)

    for l in 0:k

        res += sin(prefac * (l+1) * i) * sin(prefac * (l+1))

    end

    return chop(2 * sin(prefac) * res / (k + 2) /  sin(prefac * i))

end

#println(flatness_op(1))
#println(flatness_op(2))
#println(flatness_op(3))
#println(flatness_op(4))

function heat_ker(g::Float64,i::Int)

    res = 0.

    prefac = pi / (k + 2)

    for l in 0:k

        expo = sin(prefac / 2 * l) * sin(prefac * (l/2 +1)) / sin(prefac)^2

        res += sin(prefac * (l+1) * i) * sin(prefac * (l+1)) *
            exp(-g^2 * expo)

    end

    return chop(2 * sin(prefac) * res / (k + 2) /  sin(prefac * i))

end

#println(heat_ker(2.,1))
#println(heat_ker(2.,2))
#println(heat_ker(2.,3))
#println(heat_ker(2.,4))
#println(heat_ker(2.,5))
