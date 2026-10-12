# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# Child of driver.py: a large chained Rosenbrock problem solved with UnoSolver.jl, meant to be
# interrupted by Ctrl+C.
#    min  sum_i (1 - x_i)^2 + 100 (x_{i+1} - x_i^2)^2
#    s.t. x_i^2 + x_{i+1}^2 <= 4,  i = 1, ..., n-1

using UnoSolver

# Behave like the REPL: Ctrl+C raises InterruptException instead of exiting the process
# (non-interactive scripts default to exit_on_sigint(true)).
Base.exit_on_sigint(false)

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200_000
const m = n - 1

function objective(::Any, x)
  f = 0.0
  @inbounds for i in 1:n-1
    f += (1.0 - x[i])^2 + 100.0 * (x[i+1] - x[i]^2)^2
  end
  return f
end

function gradient(::Any, g, x)
  fill!(g, 0.0)
  @inbounds for i in 1:n-1
    r = x[i+1] - x[i]^2
    g[i] += -2.0 * (1.0 - x[i]) - 400.0 * x[i] * r
    g[i+1] += 200.0 * r
  end
  return g
end

function constraints(::Any, c, x)
  @inbounds for i in 1:m
    c[i] = x[i]^2 + x[i+1]^2
  end
  return c
end

# Jacobian: row i has entries (i, i) and (i, i+1)
const jrows = Cint[1:m; 1:m]
const jcols = Cint[1:m; 2:n]

function jacobian(::Any, jvals, x)
  @inbounds for i in 1:m
    jvals[i] = 2.0 * x[i]
    jvals[m+i] = 2.0 * x[i+1]
  end
  return jvals
end

# Lagrangian Hessian (lower triangle): diagonal (j, j), then subdiagonal (i+1, i)
const hrows = Cint[1:n; 2:n]
const hcols = Cint[1:n; 1:m]

function hessian(::Any, hvals, x, y, rho)
  # Lagrangian = rho*f(x) - y^T c(x) (lagrangian_sign = -1 below)
  fill!(hvals, 0.0)
  @inbounds for i in 1:n-1
    hvals[i] += rho * (2.0 + 1200.0 * x[i]^2 - 400.0 * x[i+1]) - 2.0 * y[i]
    hvals[i+1] += rho * 200.0 - 2.0 * y[i]
    hvals[n+i] = -400.0 * rho * x[i]
  end
  return hvals
end

const status_names = Dict(
  UnoSolver.UNO_SUCCESS => "SUCCESS",
  UnoSolver.UNO_ITERATION_LIMIT => "ITERATION_LIMIT",
  UnoSolver.UNO_TIME_LIMIT => "TIME_LIMIT",
  UnoSolver.UNO_EVALUATION_ERROR => "EVALUATION_ERROR",
  UnoSolver.UNO_ALGORITHMIC_ERROR => "ALGORITHMIC_ERROR",
  UnoSolver.UNO_USER_TERMINATION => "USER_TERMINATION",
  UnoSolver.UNO_USER_INTERRUPTION => "USER_INTERRUPTION",
)

x0 = [isodd(i) ? -1.2 : 1.0 for i in 1:n]
stats = try
  uno("NLP", true, n, m, fill(-Inf, n), fill(Inf, n), fill(-Inf, m), fill(4.0, m),
      jrows, jcols, 2m, hrows, hcols, n + m,
      objective, constraints, gradient, jacobian, hessian, nothing, nothing, nothing,
      nothing, 'L', -1, x0;
      preset="ipopt", linear_solver="MUMPS", max_iterations=100000, time_limit=1800.0)
catch e
  if e isa InterruptException
    # Julia delivered Ctrl+C as an exception in a callback: it unwound through Uno's C++ frames
    println("ESCAPED_INTERRUPT_EXCEPTION")
    flush(stdout)
    exit(3)
  end
  rethrow()
end

println("STATUS ", get(status_names, stats.optimization_status, "UNKNOWN_$(stats.optimization_status)"))
println("ITERATIONS ", stats.iter)
if all(isfinite, stats.solution) && isfinite(stats.objective)
  println("FINITE_ITERATE")
end
flush(stdout)
