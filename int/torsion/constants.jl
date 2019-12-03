const p = 8
const k = div(p, 4) - 1
const z = (k + 1)^2
const num_q_values = p   # this should suffice
const TOLERANCE = 1.0e-13;


const A = exp(2im * pi / p);

const Asqrt = exp(1im * pi / p);
