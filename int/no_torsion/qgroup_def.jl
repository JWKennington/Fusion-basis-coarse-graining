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


w0(n::Int) = chop(sqrt(q0(2*n-1)+0im))


const w_values = Complex{Float64}[w0(n) for n in -num_q_values:num_q_values];


@inline w(n::Int) = @inbounds return w_values[n+num_q_values+1];


v0(n::Int) =  chop(minus_one_pow(n-1) * w(n));

const v_values = Complex{Float64}[v0(n) for n in -num_q_values:num_q_values];

@inline v(n::Int) = @inbounds return v_values[n+num_q_values+1];


const D = sqrt( (k+2)/2 ) / sin(pi/(k+2))


function tria(a::Int, b::Int, c::Int)
    coupling_rules(a,b,c) *
    sqrt(qf(a + b - c - 1) * qf(a - b + c - 1) * qf(-a + b + c - 1) / qf(a + b + c - 2))
end

function coupling_rules(j1::Int,j2::Int,J::Int)

    res = 0

    if (j1-1) + (j2-1) >= (J-1) &&
        (j1-1) + (J-1) >= (j2-1) &&
        (j2-1) + (J-1) >= (j1-1) &&
        (j1-1) + (j2-1) + (J-1) <= p/2. - 2

        res = 1

    end

    return res

end


function sixjr(a::Int, b::Int, e::Int, d::Int, c::Int, f::Int)

    s = 0.0

    lowerBound = max(a + b + e - 3, a + c + f - 3, b + d + f - 3, d + c + e - 3)
    upperBound = min(a + b + c + d - 4, a + d + e + f - 4, b + c + e + f - 4)

    for n in lowerBound:upperBound
        s += minus_one_pow(n) * qf(n + 1) / qf(n - a - b - e + 3) /
        qf(n - a - c - f + 3) / qf(n - b - d - f + 3) / qf(n - d - c - e + 3) /
        qf(a + b + c + d - n - 4) / qf(a + d + e + f - n - 4) /
        qf(b + c + e + f - n - 4)
    end

    s *= sqrt(q(2 * e - 1) * q(2 * f - 1)) *
    minus_one_pow(a + b + c + d - 4) * tria(a, b, e) * tria(c, d, e) *
    tria(a, c, f) * tria(b, d, f)

    return s

end



function Ribbon_exp(i::Int,j::Int,a::Int,b::Int)

    #i and j are 1/2 integers

    klev = div(p,2) - 2 # This is the level k of the quantum group

    prefac = pi / (klev + 2) # This prefactor appears in all arguments

    eigval1 = sin(prefac) * sin( prefac * i * (2 * a - 1)) /
        sin(prefac * i) / sin(prefac * (2 * a -1))

    eigval2 = sin(prefac) * sin( prefac * j * (2 * b - 1)) /
            sin(prefac * j) / sin(prefac * (2 * b -1))

    res = eigval1 * eigval2

    res

end


function diff_from_id(a::Int)

    res = 0.

    for i in 1:k+1

        res += 1/v(i)^4 * (Ribbon_exp(i,i,a,a) - Ribbon_exp(i,i,1,1))

    end

    chop(real(res))

end


function Wilson(g::Float64,i::Int)

    #Define the action by comparing eigenvalues of the Wilson loop operator
    #for label l=1/2

    klev = div(p,2) - 2

    prefac = pi / (klev + 2)

    eigval = sin(prefac) * sin(2 * prefac * (2 * i - 1)) /
        sin(2 * prefac) / sin(prefac * (2 * i -1))

    res = exp(1 / g^2 * (eigval - 1))

    res

end

#println(div(p,2) - 2)

#println(Wilson(2.,1))
#println(Wilson(2.,2))
#println(Wilson(2.,3))

#Code the heat kernel action

function flatness_op(i::Int)

    res = 0.

    klev = div(p,2) - 2

    prefac = pi / (klev + 2)

    for l in 0:klev

        res += sin(prefac * (l+1) * (2*i-1)) * sin(prefac * (l+1))

    end

    return chop(2 * sin(prefac) * res / (klev + 2) /  sin(prefac * (2*i-1)))

end

#println(flatness_op(1))
#println(flatness_op(2))
#println(flatness_op(3))
#println(flatness_op(4))

function heat_ker(g::Float64,i::Int)

    res = 0.

    klev = div(p,2) - 2

    prefac = pi / (klev + 2)

    for l in 0:klev

        expo = sin(prefac / 2 * l) * sin(prefac * (l/2 +1)) / sin(prefac)^2

        res += sin(prefac * (l+1) * (2*i-1)) * sin(prefac * (l+1)) *
            exp(-g^2 * expo)

    end

    return chop(2 * sin(prefac) * res / (klev + 2) /  sin(prefac * (2*i-1)))

end

#println(heat_ker(2.,1))
#println(heat_ker(2.,2))
#println(heat_ker(2.,3))
#println(heat_ker(2.,4))
