# UnoSolver.jl

UnoSolver.jl is a wrapper for [Uno](https://github.com/cvanaret/Uno), a modern and modular solver for nonlinearly constrained optimization.

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

## Examples

Below are examples showing how to use `UnoSolver.jl` with the interfaces for `NLPModels.jl` and `MathOptInterface.jl`.

```julia
using UnoSolver, CUTEst

nlp = CUTEstModel{Float64}("HS15")
model, solver = uno(nlp, preset="filtersqp", print_solution=true, logger="INFO")

timer = UnoSolver.uno_get_cpu_time(solver)
niter = UnoSolver.uno_get_number_iterations(solver)
nsub = UnoSolver.uno_get_number_subproblem_solved_evaluations(solver)
optimization_status = UnoSolver.uno_get_optimization_status(solver)
solution_status = UnoSolver.uno_get_solution_status(solver)
solution_objective = UnoSolver.uno_get_solution_objective(solver)
solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)

primal_solution = Vector{Float64}(undef, nlp.meta.nvar)
UnoSolver.uno_get_primal_solution(solver, primal_solution)

constraint_dual_solution = Vector{Float64}(undef, nlp.meta.ncon)
UnoSolver.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

lower_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

upper_bound_dual_solution = Vector{Float64}(undef, nlp.meta.nvar)
UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)
```

```julia
using UnoSolver, JuMP

jump_model = Model(() -> UnoSolver.Optimizer(preset="filtersqp"))
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

If you encounter any issues with the interface for JuMP problems, please [open an issue](https://github.com/cvanaret/Uno/issues) so we can fix it.
As a temporary workaround, you can use [NLPModelsJuMP.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl) to wrap a JuMP model into a `MathOptNLPModel`:

```julia
using UnoSolver, NLPModelsJuMP

nlp = MathOptNLPModel(jump_model)

model = uno_model(nlp)
solver = uno_solver()
uno_set_solver_preset(solver, "filtersqp")
uno_set_solver_bool_option(solver, "print_solution", true)
uno_optimize(solver, model)
```

## Linear solvers

`UnoSolver.jl` supports a number of linear solvers. If not specified by the user, the default linear solver is picked in this order (if available): MA57, MA27, MUMPS.

### LibHSL

We highly recommend downloading the latest release of [libHSL](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL) and installing the official version of `HSL_jll.jl` into your current environment using:
```julia
import Pkg
Pkg.develop(path = "/full/path/to/HSL_jll.jl")
```

This optional dependency provides access to more reliable and powerful linear solvers. Currently, `UnoSolver.jl` supports `MA27` and `MA57`.
Pick a linear solver by setting the `linear_solver` attribute:
```julia
using JuMP, UnoSolver
import HSL_jll
model = Model(() -> UnoSolver.Optimizer(preset="ipopt"))
set_attribute(model, "linear_solver", "MA57")
```

### MUMPS

MUMPS can be used by setting the `linear_solver` attribute:
```julia
using JuMP, UnoSolver
model = Model(() -> UnoSolver.Optimizer(preset="ipopt"))
set_attribute(model, "linear_solver", "MUMPS")
```

## QP solvers

If not specified by the user, the default QP solver is BQPD.

### BQPD

BQPD can be used by setting the `QP_solver` attribute:
```julia
using JuMP, UnoSolver
model = Model(() -> UnoSolver.Optimizer(preset="filtersqp"))
set_attribute(model, "QP_solver", "BQPD")
```

## LP solvers

If not specified by the user, the default LP solver is BQPD.

### BQPD

BQPD can be used by setting the `LP_solver` attribute:
```julia
using JuMP, UnoSolver
model = Model(() -> UnoSolver.Optimizer(preset="filterslp"))
set_attribute(model, "LP_solver", "BQPD")
```

### HiGHS

HiGHS can be used by setting the `LP_solver` attribute:
```julia
using JuMP, UnoSolver
model = Model(() -> UnoSolver.Optimizer(preset="filterslp"))
set_attribute(model, "LP_solver", "HiGHS")
```

## BLAS and LAPACK demuxer

`Uno_jll.jl` is compiled with [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT), a library that can switch between BLAS and LAPACK backends at runtime, such as OpenBLAS, Intel MKL, and Apple Accelerate.
The default BLAS and LAPACK backend used in the Julia interface `UnoSolver.jl` is [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS).

### Display backends

You can check which backends are currently loaded with:

```julia
import LinearAlgebra
LinearAlgebra.BLAS.lbt_get_config()
```
If no `BLAS` or `LAPACK` library compiled with 32-bit integers (`LP64`) is available, `UnoSolver.jl` will automatically load a compatible version of `OpenBLAS`.
You can run the command again after `using UnoSolver` to verify which backend is in use.

### Sequential BLAS and LAPACK

If you have the LP64 reference versions of [BLAS and LAPACK](https://github.com/Reference-LAPACK/lapack) installed, you can switch to the sequential backends by running:

```julia
using ReferenceBLAS32_jll, LAPACK32_jll
LinearAlgebra.BLAS.lbt_forward(libblas32)
LinearAlgebra.BLAS.lbt_forward(liblapack32)
using UnoSolver
```

### MKL

If you have [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl) installed,
switch to MKL by adding `using MKL` to your code:

```julia
using MKL
using UnoSolver
```

### AppleAccelerate

If you are using macOS v13.4 or later and you have [AppleAccelerate.jl](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) installed, add `using AppleAccelerate` to your code:

```julia
using AppleAccelerate
using UnoSolver
```