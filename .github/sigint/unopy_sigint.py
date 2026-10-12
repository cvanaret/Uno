# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# Child of driver.py: a large chained Rosenbrock problem solved with unopy, meant to be interrupted
# by Ctrl+C. After the solve, checks that Python's own SIGINT handling is back in place.
#    min  sum_i (1 - x_i)^2 + 100 (x_{i+1} - x_i^2)^2
#    s.t. x_i^2 + x_{i+1}^2 <= 4,  i = 0, ..., n-2

import math
import os
import signal
import sys
import time

import numpy as np
import unopy

n = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
m = n - 1
Inf = float("inf")


def objective(x):
   return float(np.sum((1. - x[:-1])**2) + 100.*np.sum((x[1:] - x[:-1]**2)**2))


def objective_gradient(x, gradient):
   residual = x[1:] - x[:-1]**2
   gradient[:] = 0.
   gradient[:-1] = -2.*(1. - x[:-1]) - 400.*x[:-1]*residual
   gradient[1:] += 200.*residual


def constraints(x, constraint_values):
   constraint_values[:] = x[:-1]**2 + x[1:]**2


# Jacobian: row i has entries (i, i) and (i, i+1)
jacobian_row_indices = np.concatenate([np.arange(m), np.arange(m)]).tolist()
jacobian_column_indices = np.concatenate([np.arange(m), np.arange(1, n)]).tolist()


def jacobian(x, jacobian_values):
   jacobian_values[:m] = 2.*x[:-1]
   jacobian_values[m:] = 2.*x[1:]


# Lagrangian Hessian (lower triangle): diagonal (j, j), then subdiagonal (i+1, i)
hessian_row_indices = np.concatenate([np.arange(n), np.arange(1, n)]).tolist()
hessian_column_indices = np.concatenate([np.arange(n), np.arange(m)]).tolist()


def lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values):
   # Lagrangian = rho*f(x) - y^T c(x)
   diagonal = np.zeros(n)
   diagonal[:-1] += objective_multiplier*(2. + 1200.*x[:-1]**2 - 400.*x[1:])
   diagonal[1:] += objective_multiplier*200.
   diagonal[:-1] -= 2.*multipliers
   diagonal[1:] -= 2.*multipliers
   hessian_values[:n] = diagonal
   hessian_values[n:] = -400.*objective_multiplier*x[:-1]


if __name__ == "__main__":
   model = unopy.Model(unopy.PROBLEM_NONLINEAR, n, unopy.ZERO_BASED_INDEXING)
   model.set_name(f"chained_rosenbrock_{n}")
   model.set_variables_lower_bounds([-Inf]*n)
   model.set_variables_upper_bounds([Inf]*n)
   model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
   model.set_constraints(m, constraints, [-Inf]*m, [4.]*m, 2*m, jacobian_row_indices, jacobian_column_indices,
      jacobian)
   model.set_lagrangian_hessian(n + m, unopy.LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
      lagrangian_hessian)
   model.set_lagrangian_sign_convention(unopy.MULTIPLIER_NEGATIVE)
   model.set_initial_primal_iterate([-1.2 if i % 2 == 0 else 1. for i in range(n)])

   solver = unopy.UnoSolver()
   solver.set_preset("ipopt")
   solver.set_option("linear_solver", "MUMPS")
   solver.set_option("max_iterations", 100000)
   solver.set_option("time_limit", 1800.)  # safety net

   result = solver.optimize(model)
   print(f"STATUS {result.optimization_status.name}", flush=True)
   print(f"ITERATIONS {result.number_iterations}", flush=True)
   if all(math.isfinite(value) for value in result.primal_solution) and math.isfinite(result.solution_objective):
      print("FINITE_ITERATE", flush=True)

   # the solve must have restored Python's SIGINT handler: a new Ctrl+C raises KeyboardInterrupt
   try:
      os.kill(os.getpid(), signal.SIGINT)
      time.sleep(5.)
   except KeyboardInterrupt:
      print("PYTHON_HANDLER_RESTORED", flush=True)
