import numpy as np
import unopy

USE_LOGFILE = True

# Setup simple problem
n = 1
blx = np.array([-5.0])
bux = np.array([5.0])
x0 = np.array([-2.0])

m = 1
blc = np.array([1.0])
buc = np.array([float('inf')])

nnz = 1
row_indices = np.array([0], dtype=np.int32)
col_indices = np.array([0], dtype=np.int32)

model = unopy.Model(unopy.PROBLEM_NONLINEAR, n, blx, bux, unopy.ZERO_BASED_INDEXING)

def objective(x): return x[0]**2
def objective_gradient(x, grad): grad[0] = 2*x[0]
model.set_objective(unopy.MINIMIZE, objective, objective_gradient)

def constraints(x, con_val): con_val[0] = x[0]
def jacobian(x, jac_val): jac_val[0] = 1.0
model.set_constraints(m, constraints, blc, buc, nnz, row_indices, col_indices, jacobian)

model.set_initial_primal_iterate(x0)

solver = unopy.UnoSolver()
solver.set_preset("filtersqp")
solver.set_option("logger", "INFO")

if not USE_LOGFILE:
    result = solver.optimize(model)

else:
    with open("unopy_crash_flush.log", "w") as log_file:
        solver.set_logger_stream(log_file)
        result = solver.optimize(model)

print(f"Optimization finished. Status: {result.optimization_status}")