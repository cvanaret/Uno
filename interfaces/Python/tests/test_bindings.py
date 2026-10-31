# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# behavior of the Python bindings: arrays exchanged with the callbacks, argument conversion and validation, results,
# logging

import io
import numpy as np
import pytest
import unopy
from problems import Inf, hs001, hs015


@pytest.fixture
def ipopt(solver):
	solver.set_preset("ipopt")
	solver.set_option("linear_solver", "MUMPS")
	return solver


def describe(array):
	return array.dtype, array.ndim, array.flags.c_contiguous, array.flags.writeable


@pytest.mark.parametrize("number_constraints", [0, 1])
def test_callback_arrays(ipopt, number_constraints):
	"""Inputs are read-only float64 1D arrays (including empty ones), outputs are writeable"""
	seen = {}

	def record(name, *arrays):
		seen.setdefault(name, describe(arrays[0]) + tuple(len(a) for a in arrays))
		for a in arrays[1:]:
			seen.setdefault(name + " multipliers", describe(a))

	def objective(x):
		record("x", x)
		return float(x @ x)

	def objective_gradient(x, gradient):
		record("gradient", gradient)
		gradient[:] = 2.*x

	def constraints(x, values):
		record("constraints", values)
		values[0] = x[0] + x[1]

	def jacobian(x, values):
		record("jacobian", values)
		values[:] = 1.

	def lagrangian_hessian(x, objective_multiplier, multipliers, values):
		record("hessian", values)
		record("hessian input", multipliers)
		values[:] = 2.*objective_multiplier

	def notify(primals, lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers, *residuals):
		record("notify", primals, lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers)

	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	model.set_objective(unopy.MINIMIZE, objective, objective_gradient)
	if number_constraints == 1:
		model.set_constraints(1, constraints, [1.], [Inf], 2, [0, 0], [0, 1], jacobian)
	model.set_lagrangian_hessian(2, unopy.LOWER_TRIANGLE, [0, 1], [0, 1], lagrangian_hessian)
	model.set_initial_primal_iterate([3., 4.])
	ipopt.set_notify_acceptable_iterate_callback(notify)
	result = ipopt.optimize(model)
	assert result.optimization_status == unopy.OptimizationStatus.SUCCESS

	read_only, writeable = (np.float64, 1, True, False), (np.float64, 1, True, True)
	assert seen["x"] == read_only + (2,)
	assert seen["gradient"] == writeable + (2,)
	assert seen["hessian"] == writeable + (2,)
	assert seen["hessian input"] == read_only + (number_constraints,)
	# the solver callbacks receive Uno's internal iterate, which may contain slack variables (one per inequality
	# constraint with the ipopt preset), and the bound multipliers of all its variables
	number_primals = seen["notify"][4]
	assert number_primals >= 2
	assert seen["notify"] == read_only + (number_primals, number_primals, number_primals, number_constraints)
	assert seen["notify multipliers"] == read_only
	if number_constraints == 1:
		assert seen["constraints"] == writeable + (1,)
		assert seen["jacobian"] == writeable + (2,)


ARRAY_LIKES = [
	pytest.param(list, id="list"),
	pytest.param(tuple, id="tuple"),
	pytest.param(np.array, id="numpy default dtypes"),  # float64 bounds, int64 indices
	pytest.param(lambda values: np.array(values)[::-1][::-1], id="numpy view"),
]


@pytest.mark.parametrize("convert", ARRAY_LIKES)
def test_setters_accept_array_likes(ipopt, convert):
	"""All array-like inputs define the same model"""
	reference = ipopt.optimize(hs015())
	result = ipopt.optimize(hs015(convert=convert))
	assert result.number_iterations == reference.number_iterations
	assert np.array_equal(result.primal_solution, reference.primal_solution)


def test_upper_triangle_hessian(ipopt):
	lower = ipopt.optimize(hs015(triangle=unopy.LOWER_TRIANGLE))
	upper = ipopt.optimize(hs015(triangle=unopy.UPPER_TRIANGLE))
	assert upper.number_iterations == lower.number_iterations
	assert np.array_equal(upper.primal_solution, lower.primal_solution)


@pytest.mark.parametrize("setter", [
	lambda model: model.set_variables_lower_bounds([0.]),
	lambda model: model.set_variables_upper_bounds(np.zeros(3)),
	lambda model: model.set_initial_primal_iterate([0., 0., 0.]),
	lambda model: model.set_constraints(1, lambda x, c: None, [0., 1.], [1.], 2, [0, 0], [0, 1], lambda x, j: None),
	lambda model: model.set_constraints(1, lambda x, c: None, [0.], [1.], 2, [0], [0, 1], lambda x, j: None),
	lambda model: model.set_lagrangian_hessian(2, unopy.LOWER_TRIANGLE, [0, 1], np.array([0]), lambda *a: None),
], ids=["lower bounds", "upper bounds", "initial point", "constraint bounds", "Jacobian indices", "Hessian indices"])
def test_setters_reject_wrong_sizes(setter):
	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	with pytest.raises(ValueError, match="Dimension mismatch"):
		setter(model)


@pytest.mark.parametrize("indices", [np.array([0., 1.]), [0., 1.], [0.5, 1.]], ids=["numpy float64", "floats", "fractions"])
def test_indices_reject_floats(indices):
	"""Indices are never silently truncated"""
	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	with pytest.raises(TypeError):
		model.set_constraints(1, lambda x, c: None, [0.], [1.], 2, [0, 0], indices, lambda x, j: None)
	with pytest.raises(TypeError):
		model.set_lagrangian_hessian(2, unopy.LOWER_TRIANGLE, indices, [0, 1], lambda *a: None)


@pytest.mark.parametrize("value", [np.zeros((1, 2)), [[0., 0.]], "ab", [0., None]],
	ids=["2D array", "nested list", "string", "None entry"])
def test_setters_reject_non_vectors(value):
	model = unopy.Model(unopy.PROBLEM_NONLINEAR, 2, unopy.ZERO_BASED_INDEXING)
	with pytest.raises(TypeError):
		model.set_variables_lower_bounds(value)


def test_result_arrays_are_copies(ipopt):
	result = ipopt.optimize(hs001())
	solution = result.primal_solution.copy()
	result.primal_solution[:] = 0.
	assert np.array_equal(result.primal_solution, solution)


def test_logger_stream(ipopt):
	stream = io.StringIO()
	ipopt.set_option("logger", "INFO")
	ipopt.set_logger_stream(stream)
	ipopt.optimize(hs015())
	log = stream.getvalue()
	assert "Uno" in log and "Optimization status" in log
	assert len(log.splitlines()) > 10

	stream.truncate(0)
	ipopt.set_option("logger", "SILENT")
	ipopt.optimize(hs015())
	assert stream.getvalue() == ""
