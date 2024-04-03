using IntervalArithmetic, Combinatorics, Polynomials, PolynomialRoots, Serialization

prec = 1024
setprecision(prec)
K = 100
K₀ = 19
Ms = big.(collect(K:-1:0))
indices = big.([1; cumsum(Ms.+1).+1])

sups = zeros(Interval{BigFloat}, (sum(Ms.+1)))
ϵ = interval(2^(-prec//2))
println("computing sups")
for k in 0:K₀
    # println(k)
    for m=0:Ms[k+1]
        # println((k,m))
        Lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m+2*k+1,m-j)//factorial(j)) for j=0:m])
        if m>0
            lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m-1+2*k+2,m-1-j)//factorial(j)) for j=0:m-1])
        else
            lₘ = Polynomial([interval(BigFloat,0)])
        end
        P = Polynomial([interval(0), interval(2)])*lₘ+Polynomial([interval(-2*k-1), interval(1)])*Lₘ
        x = interval.(real.(PolynomialRoots.roots(mid.(P.coeffs))))
        xl = x .- ϵ
        xu = x .+ ϵ

        if all(sign.(P.(xl)).*sign.(P.(xu)) .==-1) && all(inf.(xl[1:end-1]).> sup.(xu[2:end])) && length(x) == m+1
            Z = sqrt(interval.(BigFloat, factorial(m+2*k+1)//factorial(m)//2))
            xrig = interval.(inf.(xl), sup.(xu));
            sups[indices[k+1]+m] = interval(maximum(sup.(abs.(Lₘ.(xrig).*xrig.^k .*sqrt.(xrig).*exp.(-xrig/interval(2))))))/Z
        else
            println((k, m))
        end
    end
end

serialize("suppsi", interval.(Float64, sups))