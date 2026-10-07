# Uno (R package)

[Uno](https://CRAN.R-project.org/package=Uno) is an R wrapper for Uno, available on CRAN. It provides R bindings to Uno's C API: you describe a nonlinear program with R callbacks for the objective, gradient, constraints, Jacobian and Lagrangian Hessian, and Uno solves it. It is the R analog of the [`unopy`](../Python) Python interface.

The R interface is developed and maintained by [Balasubramanian Narasimhan](https://github.com/bnaras). Sources, documentation and issue tracker live in a separate repository: [github.com/bnaras/Uno](https://github.com/bnaras/Uno) (package website: [bnaras.github.io/Uno](https://bnaras.github.io/Uno/)).

## Installation

`Uno` is a registered CRAN package and can be installed with:

```r
install.packages("Uno")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bnaras/Uno")
```

A C++17 compiler and CMake (>= 3.16) are required for source builds: the package builds Uno and the HiGHS QP/LP subproblem solver from source. The MUMPS linear solver used by the interior-point preset is reached at run time through the [rmumps](https://CRAN.R-project.org/package=rmumps) package, so no separate MUMPS installation is needed.

## Getting started

Hock--Schittkowski problem 15 (`x* = (0.5, 2)`, `f* = 306.5`), solved with the
interior-point preset. Derivatives are supplied in COO form (0-based indices);
the Hessian is the lower triangle of the Lagrangian.

```r
library(Uno)

objective   <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
gradient    <- function(x) c(400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2,
                             200 * (x[2] - x[1]^2))
constraints <- function(x) c(x[1] * x[2], x[1] + x[2]^2)
jacobian    <- function(x) c(x[2], 1, x[1], 2 * x[2])
hessian     <- function(x, sigma, lambda)
  c(sigma * (1200 * x[1]^2 - 400 * x[2] + 2),
    -400 * sigma * x[1] - lambda[1],
    200 * sigma - 2 * lambda[2])

res <- uno_solve(
  n = 2L, lb = c(-Inf, -Inf), ub = c(0.5, Inf), sense = "minimize",
  obj = objective, grad = gradient,
  m = 2L, cl = c(1, 0), cu = c(Inf, Inf), cons = constraints,
  jac_rows = c(0L, 1L, 0L, 1L), jac_cols = c(0L, 0L, 1L, 1L), jac = jacobian,
  hess_rows = c(0L, 1L, 1L), hess_cols = c(0L, 0L, 1L), hess = hessian,
  x0 = c(-2, 1), preset = "ipopt", base_indexing = 0L, verbose = FALSE,
  options = list(logger = "SILENT")
)

res$objective   # 306.5
res$primal      # 0.5 2
```

Any Uno solver option can be passed through `options` as a named list (applied
after the preset). See `vignette("Uno")` and the [package website](https://bnaras.github.io/Uno/) for the full walk-through.

If you encounter any issues with the R interface, please [open an issue](https://github.com/bnaras/Uno/issues) on the package repository.
