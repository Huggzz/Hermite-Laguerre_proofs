using IntervalArithmetic, Combinatorics, Serialization, FastGaussQuadrature, Polynomials, Base.Threads

enlarge(x::Interval{BigFloat}) = interval(BigFloat(inf(x), RoundDown), BigFloat(sup(x), RoundUp))

function get_roots(N)
    M = Threads.nthreads()
    K = N÷M+1
    ind = [1+K*m for m=0:M]
    L = Polynomial(mid.([interval((-1)^k)*interval(BigFloat, binomial(N,k)//factorial(k)) for k=0:N]))
    dL = derivative(L)
    Y, _ = gausslaguerre(Int64(N), 0.0)
    X = [big.(Y[ind[m]:min(ind[m+1]-1,N)]) for m=1:M]
    Threads.@threads for i=1:Threads.nthreads()
        x = copy(X[i])
        dx = ones(BigFloat, length(x))
        while maximum(abs.(dx))> big(2)^(-precision(BigFloat)//2)
            # println((i,Float64(log2(maximum(abs.(dx))))))
            indices = abs.(dx) .> big(2)^(-precision(BigFloat)//2)
            dx[indices] .= L.(x[indices])./dL.(x[indices])
            x[indices] .-= dx[indices]
        end
        X[i] = x
    end
    return reduce(vcat,X)
end

function laguerre(n, Y)
    M = Threads.nthreads()
    K = length(Y)÷M+1
    ind = [1+K*m for m=0:M]
    X = [Y[ind[m]:min(ind[m+1]-1,N)] for m=1:M]
    p = tuple([interval((-1)^k)*interval(BigFloat, binomial(n,k)//factorial(k)) for k=0:n])
    Threads.@threads for i=1:Threads.nthreads()
        x = X[i]
        X[i] = evalpoly.(x, p)
    end
    return reduce(vcat,X)
end

println("six-product quadrature rule:")
# all computations are bounded rigorously via Interval Arithmetic
prec = 8192
setprecision(prec)

n = big(500)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(1502)

# fast computation in arbitrary precisions of the roots
X = get_roots(N)

# rigorous enclosure of the roots of Lₙ via ϵ-inflation
ϵ = big(2)^(-prec//2)
Xl = interval.(X .- ϵ)
Xu = interval.(X .+ ϵ)
# check enclosure via intermediate value Thm and fundamental Thm of Algebra
if all(sign.(laguerre(N, Xl)).*sign.(laguerre(N, Xu)) .==-1) && all(sup.(Xl[1:end-1]).< inf.(Xu[2:end]))
    println("enclosure of roots checked")
else
    println("enclosure of roots failed")
end

Xrig = interval.(inf.(Xl), sup.(Xu));

# weights for the Gaussian-Legendre quadrature
W = Xrig ./((laguerre(N+1,Xrig)).^interval(2))./(interval(N+1)^interval(2))
W̄ = (cbrt.(sqrt.(W/interval(5))))
# quadrature is exact for quadruple product of polynomials of degree n÷2-1

# Vandermonde matrix is with respect to normalised Laguerre polynomials

setprecision(prec÷2)
V6 = zeros(Interval{BigFloat}, (N, n+1))
# change of variables
Yrig6 = Xrig ./interval(5);

Lprev = ones(Interval{BigFloat}, N) 
L = interval(BigFloat, 1) .- Yrig6
V6[:,1] = Lprev
for m=1:n
    # println(m)
    V6[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, 2*m+1) .-Yrig6).*L .- interval(BigFloat, m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end

V6 .*= enlarge.(W̄)

# quick tests
if in_interval(1349//15625, sum(V6[:,2].^interval(6)))
    println("integral ∫(1-z)^6 exp(-5z) dz = 1349/15625 checked")
else
    println("integral verification failed")
end

# save matrix

setprecision(256)

serialize("V6r", enlarge.(V6))

println("four-product quadrature rule:")
# all computations are bounded rigorously via Interval Arithmetic
prec = 4096
setprecision(prec)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2n-1
N = big(1002)
# fast computation in arbitrary precisions of the roots
X = get_roots(N)

# rigorous enclosure of the roots of Lₙ via ϵ-inflation
ϵ = big(2)^(-prec//2)
Xl = interval.(X .- ϵ)
Xu = interval.(X .+ ϵ)
# check enclosure via intermediate value Thm and fundamental Thm of Algebra
if all(sign.(laguerre(N, Xl)).*sign.(laguerre(N, Xu)) .==-1) && all(sup.(Xl[1:end-1]).< inf.(Xu[2:end]))
    println("enclosure of roots checked")
else
    println("enclosure of roots failed")
end

Xrig = interval.(inf.(Xl), sup.(Xu));

# weights for the Gaussian-Legendre quadrature
W = Xrig ./((laguerre(N+1,Xrig)).^interval(2))./(interval(N+1)^interval(2))
W̄ = (sqrt.(sqrt.(W/interval(3))))

setprecision(prec ÷ 2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V4 = zeros(Interval{BigFloat}, (N, n+1))
# change of variables
Yrig4 = Xrig ./interval(3);

Lprev = ones(Interval{BigFloat}, N) 
L = interval(BigFloat, 1) .- Yrig4
V4[:,1] = Lprev
for m=1:n
    # println(m)
    V4[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, 2*m+1) .-Yrig4).*L .- interval(BigFloat, m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end

V4 .*= enlarge.(W̄)

# quick tests
if in_interval(11//81,  sum(V4[:,2].^interval(4)))
    println("integral ∫(1-z)^4 exp(-3z) dz=11/81 checked")
else
    println("integral verification failed")
end

# save matrix

setprecision(256)

serialize("V4r", enlarge.(V4))