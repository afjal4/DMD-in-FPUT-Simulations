using SciMLOperators, DifferentialEquations, LinearAlgebra, LinearSolve, SciMLOperators

function shift_left!(w, v, u, N, t)
    @inbounds for i in 1:N-1
        w[i] = v[i+1]
    end
    w[N] = 0.0
    nothing
end

function shift_right!(y, x, N, t)
    @inbounds for i in 1:N-1
        y[i+1] = x[i]
    end
    y[1] = 0.0
end

function shift_right!(v, u, N, t)
    w = zeros(N)
    shift_left!(w, v, u, N, t)
    w
end

N = 10000
x = rand(N)

u = zeros(N)
p = N
t= 0.0

l_shift_op = FunctionOperator(shift_left!, zeros(N), zeros(N); u,p,t) #!!agh why does this not work?

b = l_shift_op * x

#=
function jacobian_vector_prod!(Jv, v; α=1.0)
    N = length(u) ÷ 2
    Jv = similar(u)

    # dF/dq = 1
    Jv[1:N] .= v[N+1:2N]


    # dF/dp = f'(q) = A(1 + Bαq) + BAαq, f(q) = (Aq).(1 + αBq) 
    function A(v)
        return circshift(v, -1) - 2 * v + circshift(v, 1)
    end

    function B(v)
        return circshift(v, -1) - circshift(v, 1)
    end
    
    q = v[1:N]
    Jv[N+1:2N] .= A(1 .+ B(α * q)) + α * B(A(q))

    return Jv
end
=#