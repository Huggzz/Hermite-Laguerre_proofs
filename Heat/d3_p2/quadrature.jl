using IntervalArithmetic, Combinatorics, Serialization, FastGaussQuadrature, Polynomials, Base.Threads

function enlarge(x::Interval{BigFloat})
    isguaranteed(x) || error("interval is not guaranteed")
    return interval(BigFloat(inf(x), RoundDown), BigFloat(sup(x), RoundUp))
end

function get_roots(N)
    M = Threads.nthreads()
    K = N÷M+1
    ind = [1+K*m for m=0:M]
    L = Polynomial(mid.([interval((-1)^k)*interval(BigFloat, doublefactorial(2*N+1)//doublefactorial(2*k+1)//factorial(k)//factorial(N-k)//2^(N-k)) for k=0:N]))
    dL = derivative(L)
    Y = ((FastGaussQuadrature.hermite_initialguess(2*Int64(N)+1)).^2)[2:end]
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
    p =  tuple([interval((-1)^k)*interval(BigFloat, doublefactorial(2*n+1)//doublefactorial(2*k+1)//factorial(k)//factorial(n-k)//2^(n-k)) for k=0:n])
    Threads.@threads for i=1:Threads.nthreads()
        x = X[i]
        X[i] = evalpoly.(x, p)
    end
    return reduce(vcat,X)
end

prec = 4096
setprecision(prec)
n = big(500)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(1002)

# fast computation in arbitrary precisions of the roots of L_N
X = get_roots(N)

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
W = sqrt(interval(BigFloat,π))*interval(BigFloat, doublefactorial(2*N+1)//factorial(N)//2^(N+1)).*Xrig ./((laguerre(N+1, Xrig)).^interval(2))./interval((N+1)^2);

setprecision(prec÷2)
# Vandermonde matrix is with respect to normalised Laguerre polynomials
V4 = zeros(Interval{BigFloat}, (N, n +1));
V3 = zeros(Interval{BigFloat}, (N, n + 1));
# change of variables
Yrig4 = enlarge.(Xrig ./interval(3));
Yrig3 = enlarge.(Xrig ./interval(2));

# computing pseudo-vandermonde matrices
println("computing Vandermonde matrices")
Zs = [sqrt(sqrt(interval(BigFloat,π))*interval(BigFloat, doublefactorial(2*j+1)//2^(j+1)//factorial(j))) for j=0:n]
α = 1//2
Lprev = ones(Interval{BigFloat}, length(Yrig3)) 
L = Lprev .+ interval(BigFloat,α).- Yrig3
V3[:, 1] = Lprev/Zs[1]
for k=1:n
    V3[:,k+1] = L/Zs[k+1]
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig3/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev)
    Lprev[:] = Lhold
end
Lprev = ones(Interval{BigFloat}, length(Yrig4)) 
L = Lprev .+ interval(BigFloat,α).- Yrig4
V4[:, 1] = Lprev/Zs[1]
for k=1:n
    V4[:,k+1] = L/Zs[k+1]
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig4/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev)
    Lprev[:] = Lhold
end

# regularized Vandermonde matrix
V̄3 = enlarge.(cbrt.(W*interval(BigFloat,2)^interval(-3//2))).*V3;
V̄4 = enlarge.(sqrt.(sqrt.(W*interval(BigFloat,3)^interval(-3//2)))).*V4;

# quick tests
if in_interval(1969//2048 , sum(V̄3[:,3].^interval(3))*(sqrt(interval(30))*interval(π)^interval(1//4))) && in_interval(974//2187, sum(V̄4[:,3].^interval(4))*sqrt(interval(3)*interval(π)))
    println("checks on some basic integrals passed")
else
    println("integral checks failed")
end

# save matrices
setprecision(256)
serialize("V3r", enlarge.(V̄3))
serialize("V4r", enlarge.(V̄4))