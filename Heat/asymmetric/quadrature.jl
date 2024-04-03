using IntervalArithmetic, Combinatorics, Serialization, FastGaussQuadrature, Polynomials, Base.Threads, Random

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

println("four-product quadrature rule")
prec = 2048
setprecision(prec)
# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1

K = big(150)
N = big(2*K+2)

X = get_roots(N)

ϵ = big(2)^(-prec//2)
Xl = interval.(BigFloat, X .- ϵ)
Xu = interval.(BigFloat, X .+ ϵ)

# check enclosure via intermediate value Thm and fundamental Thm of Algebra
if all(sign.(laguerre(N,Xl)).*sign.(laguerre(N,Xu)) .==-1) && all(sup.(Xl[1:end-1]).< inf.(Xu[2:end])) && length(Xl) == N
    println("enclosure of roots checked")
else
    println("enclosure of roots failed")
end

Xrig = interval.(inf.(Xl), sup.(Xu));

# weights for the Gaussian-Laguerre quadrature
W = Xrig ./((laguerre(N+1, Xrig)).^2)./(interval(N+1)^2);
Ms = big.(collect(K:-1:0))
indices = big.([1; cumsum(Ms.+1).+1])
# Vandermonde matrix is with respect to normalised Laguerre polynomials
setprecision(prec÷2)
W̄ = enlarge.(sqrt.(sqrt.(W/interval(6))))
V4 = zeros(Interval{BigFloat}, (N, sum(Ms.+1)))
# change of variables
Yrig4 = enlarge.(Xrig ./interval(3));

ind = shuffle(collect(0:K))
# computing pseudo-vandermonde matrix
println("computing Vandermonde matrix")
Threads.@threads for k in ind
    # println(k)
    p = Yrig4 .^interval(k) .*sqrt.(Yrig4)
    Lprev = ones(Interval{BigFloat}, N) 
    L = interval(BigFloat, 2*k+2) .- Yrig4
    V4[:,indices[k+1]] =p .* W̄ .* Lprev/sqrt(interval(BigFloat, factorial(2*k+1)//2))
    for m=1:Ms[k+1]
        Z = sqrt(interval(BigFloat, factorial(m+2*k+1)//factorial(m)//2))
        V4[:,indices[k+1]+m] = p .*W̄ .* L/Z
        # println((k,m))
        L, Lprev = ((interval(BigFloat, 2*(m+k)+2) .-Yrig4).*L .- interval(BigFloat, m+2*k+1)*Lprev)/interval(BigFloat, m+1), L
    end
end

# quick tests
if in_interval(220//19683//2, sum((V4[:,indices[2]+1]/sqrt(interval(BigFloat,2))).^4))
    println("integral verification succeeded")
else
    println("integral verification failed")
end

# save matrix
serialize("V4r", interval.(Float64, V4))

println("six-product quadrature rule:")

prec = 2048
setprecision(prec)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(3*K+2)

X = get_roots(N)

ϵ = big(2)^(-prec//2)
Xl = interval.(BigFloat, X .- ϵ)
Xu = interval.(BigFloat, X .+ ϵ)

# check enclosure via intermediate value Thm and fundamental Thm of Algebra
if all(sign.(laguerre(N,Xl)).*sign.(laguerre(N,Xu)) .==-1) && all(sup.(Xl[1:end-1]).< inf.(Xu[2:end])) && length(Xl) == N
    println("enclosure of roots checked")
else
    println("enclosure of roots failed")
end

Xrig = interval.(inf.(Xl), sup.(Xu));

# weights for the Gaussian-Laguerre quadrature
W = Xrig ./((laguerre(N+1, Xrig)).^2)./(interval(N+1)^2);

Ms = big.(collect(K:-1:0))
indices = big.([1; cumsum(Ms.+1).+1])
setprecision(prec÷2)
W̄ = enlarge.(cbrt.(sqrt.(W/interval(10))))
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V6 = zeros(Interval{BigFloat}, (N, sum(Ms.+1)))
# change of variables
Yrig6 = enlarge.(Xrig ./interval(5));

ind = shuffle(collect(0:K))
# computing pseudo-vandermonde matrix
println("computing Vandermonde matrix")
Threads.@threads for k in ind
    # println(k)
    p = Yrig6.^interval(k) .*sqrt.(Yrig6)
    Lprev = ones(Interval{BigFloat}, N) 
    L = interval(BigFloat, 2*k+2) .- Yrig6
    V6[:,indices[k+1]] = p .* W̄ .* Lprev/sqrt(interval(BigFloat, factorial(2*k+1)//2))
    for m=1:Ms[k+1]
        Z = sqrt(interval(BigFloat, factorial(m+2*k+1)//factorial(m)//2))
        V6[:,indices[k+1]+m] = p .* W̄ .* L/Z
        # println((k,m))
        L, Lprev = ((interval(BigFloat, 2*(m+k)+2) .-Yrig6).*L .- interval(BigFloat, m+2*k+1)*Lprev)/interval(BigFloat, m+1), L
    end
end

# quick tests
if in_interval(532308//1220703125//2, sum((V6[:,indices[2]+1]/sqrt(interval(BigFloat,2))).^6))
    println("integral verification succeeded")
else
    println("integral verification failed")
end

serialize("V6r", interval.(Float64, V6))