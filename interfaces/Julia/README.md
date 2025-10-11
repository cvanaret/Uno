# Uno.jl

Uno.jl is a wrapper for [Uno](https://github.com/cvanaret/Uno), a modern and modular solver for nonlinearly constrained optimization.

The package has three components:

* a thin wrapper around the complete C API,
* an interface to [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) for solving any optimization problem following the API, such as [CUTEst](https://github.com/JuliaSmoothOptimizers/CUTEst.jl) problems,
* an interface to [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl) for handling [JuMP](https://github.com/jump-dev/JuMP.jl) models.

## Affiliation

This Julia interface is developed and maintained by [Alexis Montoison](https://github.com/amontoison) and [Charlie Vanaret](https://github.com/cvanaret).

## Installation

`Uno.jl` is not yet a registered Julia package, but it can still be installed and tested through the Julia package manager.

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/cvanaret/Uno", subdir="interfaces/Julia")
julia> Pkg.test("Uno")
```

We plan to register `Uno.jl` in the near future (at the latest by JuMP-dev 2025).

## Examples

Below are examples showing how to use `Uno.jl` with the interfaces for `NLPModels.jl` and `MathOptInterface.jl`.

```julia
using Uno, CUTEst

nlp = CUTEstModel{Float64}("HS15")
model = uno_model(nlp)
solver = uno_solver("filtersqp", print_solution="yes")
uno_optimize(solver, model)

timer = Uno.uno_get_cpu_time(solver)
niter = Uno.uno_get_number_iterations(solver)
nsub = Uno.uno_get_number_subproblem_solved_evaluations(solver)
optimization_status = Uno.uno_get_optimization_status(solver)
solution_status = Uno.uno_get_solution_status(solver)
solution_objective = Uno.uno_get_solution_objective(solver)
solution_primal_feasibility = Uno.uno_get_solution_primal_feasibility(solver)
solution_dual_feasibility = Uno.uno_get_solution_dual_feasibility(solver)
solution_complementarity = Uno.uno_get_solution_complementarity(solver)

primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_primal_solution(solver, primal_solution)

constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
Uno.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
Uno.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
```

```julia
using Uno, JuMP

jump_model = Model(() -> Uno.Optimizer(preset="filtersqp"))
x0 = [-2, 1]
uvar = [0.5, Inf]
@variable(jump_model, x[i = 1:2] ≤ uvar[i], start = x0[i])
@objective(jump_model, Min, 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)
@constraint(jump_model, x[1] * x[2] - 1 ≥ 0)
@constraint(jump_model, x[1] + x[2]^2 ≥ 0)

optimize!(jump_model)

termination_status(jump_model)  # solver termination status
objective_value(jump_model)     # objective value
value.(x)                       # primal solution
```

If you encounter any issues with the interface for JuMP problems, please [open an issue](https://github.com/JuliaSmoothOptimizers/Uno.jl/issues) so we can fix it.
As a temporary workaround, you can use [NLPModelsJuMP.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl) to wrap a JuMP model into a `MathOptNLPModel`:

```julia
using NLPModelsJuMP

nlp = MathOptNLPModel(jump_model)

model = uno_model(nlp)
solver = uno_solver("filtersqp", print_solution="yes")
uno_optimize(solver, model)
```

## LibHSL

We highly recommend downloading the latest release of [libHSL](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL) and installing the official version of `HSL_jll.jl`.
This optional dependency provides access to more reliable and powerful linear solvers in `Uno.jl`, such as `MA27` and `MA57`.

## BLAS and LAPACK demuxer

`Uno_jll.jl` is compiled with [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT), a library that can switch between BLAS and LAPACK backends at runtime, such as OpenBLAS, Intel MKL, and Apple Accelerate.
The default BLAS and LAPACK backend used in the Julia interface `Uno.jl` is [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS).

### Display backends

You can check which backends are currently loaded with:

```julia
import LinearAlgebra
LinearAlgebra.BLAS.lbt_get_config()
```
If no `BLAS` or `LAPACK` library compiled with 32-bit integers (`LP64`) is available, `Uno.jl` will automatically load a compatible version of `OpenBLAS`.
You can run the command again after `using Uno` to verify which backend is in use.

### Sequential BLAS and LAPACK

If you have the LP64 reference versions of [BLAS and LAPACK](https://github.com/Reference-LAPACK/lapack) installed, you can switch to the sequential backends by running:

```julia
using ReferenceBLAS32_jll, LAPACK32_jll
LinearAlgebra.BLAS.lbt_forward(libblas32)
LinearAlgebra.BLAS.lbt_forward(liblapack32)
using Uno
```

### MKL

If you have [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl) installed,
switch to MKL by adding `using MKL` to your code:

```julia
using MKL
using Uno
```

### AppleAccelerate

If you are using macOS v13.4 or later and you have [AppleAccelerate.jl](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) installed, add `using AppleAccelerate` to your code:

```julia
using AppleAccelerate
using Uno
```
