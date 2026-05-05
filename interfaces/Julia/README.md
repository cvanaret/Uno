# UnoSolver.jl

UnoSolver.jl is a wrapper for [Uno](https://github.com/cvanaret/Uno), a modern and modular solver for nonlinearly constrained optimization.
Uno unifies Lagrange-Newton (**SQP** and **interior-point**) methods by breaking them down into a set of common building blocks such as constraint reformulation, step computation, and globalization.

The package has three components:

* a thin wrapper around the complete C API,
* an interface to [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) for solving any optimization problem following the API, such as [CUTEst](https://github.com/JuliaSmoothOptimizers/CUTEst.jl) problems,
* an interface to [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl) for handling [JuMP](https://github.com/jump-dev/JuMP.jl) models.

## Affiliation

This Julia interface is developed and maintained by [Alexis Montoison](https://github.com/amontoison) and [Charlie Vanaret](https://github.com/cvanaret).

## Installation

`UnoSolver.jl` is a registered Julia package, it can be installed and tested through the Julia package manager.

```julia
julia> using Pkg
julia> Pkg.add("UnoSolver")
julia> Pkg.test("UnoSolver")
```

## Getting started

To get started with Uno, check out the [official documentation](https://unosolver.readthedocs.io/en/latest/interfaces/julia/).

## How to cite Uno

Our Uno paper was accepted in the Mathematical Programming Computation journal on Feb 22, 2026.
Our preprint is available on [ResearchGate](https://www.researchgate.net/publication/397446552_Implementing_a_unified_solver_for_nonlinearly_constrained_optimization) and can be cited with the following BibTeX entry:

```
@unpublished{VanaretLeyffer2026,
  author = {Vanaret, Charlie and Leyffer, Sven},
  title = {Implementing a unified solver for nonlinearly constrained optimization},
  year = {2026},
  note = {Accepted to Mathematical Programming Computation on Feb 22, 2026}
}
```