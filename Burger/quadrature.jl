using IntervalArithmetic, Combinatorics, Serialization, FastGaussQuadrature, Polynomials, Base.Threads

function enlarge(x::Interval{BigFloat})
    isguaranteed(x) || error("interval is not guaranteed")
    return interval(BigFloat(inf(x), RoundDown), BigFloat(sup(x), RoundUp))
end


# redefine double factorial
function doublefac(n)
    if n==-1
        return 1
    else
        return doublefactorial(n)
    end
end

println("six-product quadrature rule:")
function get_roots(N)
    M = Threads.nthreads()
    K = N÷M+1
    ind = [1+K*m for m=0:M]
    L = Polynomial(mid.([interval((-1)^k)*interval(BigFloat, doublefactorial(2*N+1)//doublefactorial(2*k+1)//factorial(k)//factorial(N-k)//2^(N-k)) for k=0:N]))
    dL = derivative(L)
    # Y = ((FastGaussQuadrature.hermite_initialguess(2*Int64(N)+1)).^2)[2:end]
    Y, _ = gausslaguerre(Int64(N), 0.5)
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

prec = 16384
setprecision(prec)
n = big(1500)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(3*n+3)

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
W̄ = sqrt.(cbrt.(W*interval(BigFloat,5)^interval(-3//2)))

setprecision(prec÷4)
V6 = zeros(Interval{BigFloat}, (N, n +1));
DV6 = zeros(Interval{BigFloat}, (N, n +1));
# change of variables
Yrig6 = enlarge.(Xrig ./interval(5));

Zs = sqrt.(sqrt(interval(BigFloat,π)).*[interval(BigFloat, doublefac(2*j-1)//2^(j)//factorial(j)) for j=0:n])

α = -1//2
β = 1//2
lprev = ones(Interval{BigFloat}, length(Yrig6)) 
l = lprev .+ interval(BigFloat,β).- Yrig6
Lprev = ones(Interval{BigFloat}, length(Yrig6)) 
L = Lprev .+ interval(BigFloat,α).- Yrig6
println("computing Vandermonde matrix")
V6[:,1] = Lprev/Zs[1]
DV6[:,1] = -(Lprev)/Zs[1]
for k=1:n
    V6[:,k+1] = L/Zs[k+1]
    DV6[:,k+1] = -(lprev+L)/Zs[k+1]
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig6/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev)
    Lprev[:] = Lhold
    lhold = copy(l)
    l[:] = ((interval(BigFloat, (2*k+1+β)//(k+1)) .-Yrig6/interval(BigFloat, k+1)).*l .- interval(BigFloat, (k+β)//(k+1))*lprev)
    lprev[:] = lhold
end

V6 .*= enlarge.(W̄)
DV6 .*= enlarge.(W̄)

# quick tests
if in_interval(669//31250, sum(V6[:,2].^4 .*DV6[:,2].^2)*interval(BigFloat,π)*sqrt(interval(BigFloat, 5)))
    println("integral sqrt(5)*π∫ψ̂₁(r)⁴*ψ̂₁'(r)² exp(-r²/4) = 669/31250 checked")
else
    println("integral verification failed")
end

serialize("V6r", interval.(Float64, V6))
serialize("DV6r", interval.(Float64, DV6))

println("four-product quadrature rule")
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

prec = 16384
setprecision(prec)
# Gaussian-Laguerre quadrature exact for polynomials of degree 2N-1
N = big(2*n+3)

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

W = Xrig ./((laguerre(N+1,Xrig)).^interval(2))./(interval(N+1)^interval(2))
W̄ = sqrt.(sqrt.(W/interval(3)))

setprecision(prec÷4)
V4 = zeros(Interval{BigFloat}, (N, n +1));
DV4 = zeros(Interval{BigFloat}, (N, n +1));
Yrig4 = enlarge.(Xrig ./interval(3));

Zs = sqrt.(sqrt(interval(BigFloat,π)).*[interval(BigFloat, doublefac(2*j-1)//2^(j)//factorial(j)) for j=0:n])

α = -1//2
β = 1//2
lprev = ones(Interval{BigFloat}, length(Yrig4)) 
l = lprev .+ interval(BigFloat,β).- Yrig4
Lprev = ones(Interval{BigFloat}, length(Yrig4)) 
L = Lprev .+ interval(BigFloat,α).- Yrig4
println("computing Vandermonde matrix")
V4[:,1] = Lprev/Zs[1]
DV4[:,1] = -(Lprev)/Zs[1]
for k=1:n
    V4[:,k+1] = L/Zs[k+1]
    DV4[:,k+1] = -(lprev+L)/Zs[k+1]
    Lhold = copy(L)
    L[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig4/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev)
    Lprev[:] = Lhold
    lhold = copy(l)
    l[:] = ((interval(BigFloat, (2*k+1+β)//(k+1)) .-Yrig4/interval(BigFloat, k+1)).*l .- interval(BigFloat, (k+β)//(k+1))*lprev)
    lprev[:] = lhold
end

V4 .*= enlarge.(W̄)
DV4 .*= enlarge.(W̄)

if in_interval(-29//324, sum(V4[:,2].^3 .*DV4[:,2])*interval(BigFloat,π))
    println("integral π∫ψ̂₁(r)³*ψ̂₁'(r) exp(-r²/4) = 29/324 checked")
else
    println("integral verification failed")
end

serialize("V4r", interval.(Float64, V4))
serialize("DV4r", interval.(Float64, DV4))