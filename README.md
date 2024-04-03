This repository contains the Julia code...

Each subfolder corresponds to an equaiton studided in the paper and contains a numerical solution `ubar`, a julia script computing the necessary rigorous quadratures `quadrature.jl` which can be executed from the subfolder in the julia REPL wiht the command

```julia
    include("quadrature.jl")
```

The proof can then be executed step by step in the jupyter notebook `proof.ipynb`.

The computations in all of the quadratures and in some of the proofs have a multithreaded implementation. The Julia REPL can be started e.g. with eight threads with the command line

    julia --threads=8

This implementation crucially relies on the package [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) for rigorous float-point arithmetic.