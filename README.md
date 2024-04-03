This repository contains the Julia code...

You will need to install the relevant packages from the julia REPL with the command

```julia
    using Pkg; Pkg.add("Combinatorics, FastGaussQuadrature, IntervalArithmetic, LaTeXStrings, Plots, Polynomials, PolynomialRoots, Random, Serialization")
```

Each subfolder corresponds to an equation studied in the paper and contains a numerical solution `ubar`, a julia script computing the necessary rigorous quadratures `quadrature.jl` which can be executed from the relevant subfolder in the julia REPL with the command

```julia
    include("quadrature.jl")
```

The proof can then be executed step by step in the jupyter notebook `proof.ipynb`.

The computations in all of the quadratures and in some of the proofs have a multithreaded implementation. The Julia REPL can be started e.g. with eight threads with the command line

    julia --threads=8


This implementation crucially relies on the package [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) for rigorous float-point arithmetic.
