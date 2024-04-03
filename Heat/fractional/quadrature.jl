using IntervalArithmetic, Combinatorics, Serialization, FastGaussQuadrature, Polynomials, Base.Threads
function enlarge(x::Interval{BigFloat})
    isguaranteed(x) || error("interval is not guaranteed")
    return interval(BigFloat(inf(x), RoundDown), BigFloat(sup(x), RoundUp))
end

function get_roots(N)
    M = Threads.nthreads()
    K = N÷M+1
    ind = [1+K*m for m=0:M]
    L = Polynomial(mid.([interval((-1)^k)*interval(BigFloat, binomial(N,k)//factorial(k)) for k=0:N]))
    dL = derivative(L)
    Y, _ = gausslaguerre(Int64(N), 0.0)
    X = [big.(Y[ind[m]:min(ind[m+1]-1,N)]) for m=1:M]
    Threads.@threads for i=1:Threads.nthreads()
    # for i =1:1
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
    X = [Y[ind[m]:min(ind[m+1]-1,length(Y))] for m=1:M]
    p = tuple([interval((-1)^k)*interval(BigFloat, binomial(n,k)//factorial(k)) for k=0:n])
    Threads.@threads for i=1:Threads.nthreads()
        x = X[i]
        X[i] = evalpoly.(x, p)
    end
    return reduce(vcat,X)
end

println("conversion quadrature rule:")
prec = 4096
setprecision(prec)
n = big(400)
# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(1201)
# fast computation in arbitrary precisions of the roots of L_N
X = get_roots(N);

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

# weights for the Gaussian-Laguerre quadrature
W = Xrig ./((laguerre(N+1,Xrig)).^2)./(interval(BigFloat, N+1)^2);
W̄ = cbrt.(sqrt.(W))

setprecision(prec÷2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V3 = zeros(Interval{BigFloat}, (N, n +1));
# change of variables
Yrig3 = enlarge.(Xrig ./interval(3));
# computing pseudo-vandermonde matrices
println("computing Vandermonde matrix")
Lprev = ones(Interval{BigFloat}, length(Yrig3)) 
L = Lprev .- Yrig3
V3[:, 1] = Lprev
for m=1:n
    V3[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat,2*m+1) .-Yrig3).*L .- interval(BigFloat,m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end

# regularized Vandermonde matrix for triple integral product
V3 .*= enlarge.(W̄);

setprecision(prec)
n = 1200

W̄ = sqrt.(W)

setprecision(prec÷2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V1 = zeros(Interval{BigFloat}, (N, n+1));
# change of variables
Yrig = enlarge.(Xrig)
# computing pseudo-vandermonde matrices
println("computing Vandermonde matrix")
Lprev = ones(Interval{BigFloat}, length(Yrig)) 
L = Lprev .- Yrig
V1[:, 1] = Lprev
for m=1:n
    V1[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat,2*m+1) .-Yrig).*L .- interval(BigFloat,m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end

# regularized Vandermonde matrix
V1 .*= enlarge.(W̄);

# save matrices
setprecision(256)
serialize("V3r", enlarge.(V3))
serialize("V1r", enlarge.(V1))
V1 = 0
V3 = 0

println("8/3-product quadrature rule:")
prec = 8192
setprecision(prec)
n = big(1200)
# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(1601)
# fast computation in arbitrary precisions of the roots of L_N
X = get_roots(N);

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

# weights for the Gaussian-Laguerredre quadrature
W = Xrig ./((laguerre(N+1,Xrig)).^2)./(interval(BigFloat, N+1)^2);
W̄ = sqrt.(sqrt.(sqrt.((W)*interval(3//5)))).^3

setprecision(prec÷2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V8 = zeros(Interval{BigFloat}, (N, n +1));
# change of variables
Yrig8 = enlarge.(Xrig .*interval(3//5));
# computing pseudo-vandermonde matrices
println("computing Vandermonde matrix")
Lprev = ones(Interval{BigFloat}, length(Yrig8)) 
L = Lprev .- Yrig8
V8[:, 1] = Lprev
for m=1:n
    V8[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat,2*m+1) .-Yrig8).*L .- interval(BigFloat,m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end


V8 .*= enlarge.(W̄);

# save matrices
setprecision(256)
serialize("V8r", enlarge.(V8))
V8 = 0

println("10/3-product quadrature rule:")
prec = 8192
setprecision(prec)
n = big(1200)
# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(2001)
# fast computation in arbitrary precisions of the roots of L_N
X = get_roots(N);

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

# weights for the Gaussian-Laguerre quadrature
W = Xrig ./((laguerre(N+1,Xrig)).^2)./(interval(BigFloat, N+1)^2);
W̄ = sqrt.((W)*interval(3//7)).^interval(3//5)

setprecision(prec÷2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V10 = zeros(Interval{BigFloat}, (N, n +1));
# change of variables
Yrig10 = enlarge.(Xrig .*interval(3//7));
# computing pseudo-vandermonde matrices
println("computing Vandermonde matrix")
Lprev = ones(Interval{BigFloat}, length(Yrig10)) 
L = Lprev .- Yrig10
V10[:, 1] = Lprev
for m=1:n
    V10[:,m+1] = L
    Lhold = copy(L)
    L[:] = ((interval(BigFloat,2*m+1) .-Yrig10).*L .- interval(BigFloat,m)*Lprev)/interval(BigFloat, m+1)
    Lprev[:] = Lhold
end

# regularized Vandermonde matrix
V10 .*= enlarge.(W̄);

# save matrices
setprecision(256)
serialize("V10r", enlarge.(V10))