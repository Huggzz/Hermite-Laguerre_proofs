using IntervalArithmetic, Combinatorics, PolynomialRoots, Serialization

function laguerre(m, x)
    α = 1//2
    L1prev = ones(Interval{BigFloat}, length(x)) 
    L1 = L1prev .+ interval(BigFloat,α).- x 
    for k=1:m-1 
        L1, L1prev = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-x/interval(BigFloat, k+1)).*L1 .- interval(BigFloat, (k+α)//(k+1))*L1prev), L1
    end
    p = tuple([interval((-1)^k)*interval(BigFloat, doublefactorial(2*m+1)//doublefactorial(2*k+1)//factorial(k)//factorial(m-k)//2^(m-k)) for k=0:m])
    L2 = evalpoly.(x, p)
    L3 = zeros(Interval{BigFloat}, length(x))
    for k=1:length(x)
        if diam(L1[k])<diam(L2[k])
            L3[k] = L1[k]
        else
            L3[k] = L2[k]
        end
    end
    return L3
end

# all computations are bounded rigorously via Interval Arithmetic
prec = 4096
setprecision(prec)
n = big(500)

# Gaussian-Laguerre quadrature exact for polynomials of degree 2n-1
N = big(1002)
coeffs = [interval((-1)^k)*interval(BigFloat, doublefactorial(2*N+1)//doublefactorial(2*k+1)//factorial(k)//factorial(N-k)//2^(N-k)) for k=0:N]
# fast computation in arbitrary precisions of the roots of Lₙ
X = real.(PolynomialRoots.roots(mid.(coeffs)));

# rigorous enclosure of the roots of Lₙ via ϵ-inflation
ϵ = big(2)^(-prec//2)
Xl = interval.(X .- ϵ)
Xu = interval.(X .+ ϵ)

# check enclosure via intermediate value Thm and fundamental Thm of Algebra
if all(sign.(laguerre(N, Xl)).*sign.(laguerre(N, Xu)) .==-1) && all(inf.(Xl[1:end-1]).> sup.(Xu[2:end]))
    println("enclosure of roots checked")
else
    println("enclosure of roots failed")
end

Xrig = interval.(inf.(Xl), sup.(Xu));

# computation of Lₘ where m = n+1
# weights for the Gaussian-Legendre quadrature
W = sqrt(interval(BigFloat,π))*interval(BigFloat, doublefactorial(2*N+1)//factorial(N)//2^(N+1)).*Xrig ./((laguerre(N+1, Xrig)).^interval(2))./interval((N+1)^2);

# quadrature is exact for quadruple product of polynomials of degree n÷2-1

# Vandermonde matrix is with respect to normalised Laguerre polynomials
V4 = zeros(Interval{BigFloat}, (N, n +1));
V3 = zeros(Interval{BigFloat}, (N, n + 1));
# change of variables
Yrig4 = Xrig ./interval(3);
Yrig3 = Xrig ./interval(2);

# computing pseudo-vandermonde matrices
println("computing Vandermonde matrices")
Zs = [sqrt(sqrt(interval(π))*interval(BigFloat, doublefactorial(2*j+1)//2^(j+1)//factorial(j))) for j=0:n]
α = 1//2
Lprev = ones(Interval{BigFloat}, length(Yrig3)) 
L = Lprev .+ interval(BigFloat,α).- Yrig3
V3[:, 1] = Lprev/Zs[1]
for k=1:n
    V3[:,k+1] = L/Zs[k+1] 
    L[:], Lprev[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig3/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev), L
end
Lprev = ones(Interval{BigFloat}, length(Yrig4)) 
L = Lprev .+ interval(BigFloat,α).- Yrig4
V4[:, 1] = Lprev/Zs[1]
for k=1:n
    V4[:,k+1] = L/Zs[k+1] 
    L[:], Lprev[:] = ((interval(BigFloat, (2*k+1+α)//(k+1)) .-Yrig4/interval(BigFloat, k+1)).*L .- interval(BigFloat, (k+α)//(k+1))*Lprev), L
end

# regularized Vandermonde matrix for triple integral product
V̄3 = (cbrt.(W*interval(2)^interval(-3//2))).*V3;
V̄4 = (sqrt.(sqrt.(W*interval(3)^interval(-3//2)))).*V4;

# quick tests
if in_interval(1969//2048 , sum(V̄3[:,3].^interval(3))*(sqrt(interval(30))*interval(π)^interval(1//4))) && in_interval(974//2187, sum(V̄4[:,3].^interval(4))*sqrt(interval(3)*interval(π)))
    println("checks on some basic integrals passed")
else
    println("integral checks failed")
end

# save matrices
setprecision(256)
V̄3 *= interval(1)
V̄4 *= interval(1)
serialize("V3r", V̄3)
serialize("V4r", V̄4)