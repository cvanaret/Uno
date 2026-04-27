# hs71
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1, 5, 5, 1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

function c_objective_hs71(userdata::Nothing, x::Vector{Float64})
  f = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
  return f
end

function c_constraints_hs71(userdata::Nothing, c::Vector{Float64}, x::Vector{Float64})
  c[1] = x[1] * x[2] * x[3] * x[4]
  c[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
  return c
end

function c_objective_gradient_hs71(userdata::Nothing, g::Vector{Float64}, x::Vector{Float64})
  g[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
  g[2] = x[1] * x[4]
  g[3] = x[1] * x[4] + 1
  g[4] = x[1] * (x[1] + x[2] + x[3])
  return g
end

function c_jacobian_hs71(userdata::Nothing, jvals::Vector{Float64}, x::Vector{Float64})
  # Constraint (row) 1
  jvals[1] = x[2] * x[3] * x[4]  # 1,1
  jvals[2] = x[1] * x[3] * x[4]  # 1,2
  jvals[3] = x[1] * x[2] * x[4]  # 1,3
  jvals[4] = x[1] * x[2] * x[3]  # 1,4
  # Constraint (row) 2
  jvals[5] = 2 * x[1]  # 2,1
  jvals[6] = 2 * x[2]  # 2,2
  jvals[7] = 2 * x[3]  # 2,3
  jvals[8] = 2 * x[4]  # 2,4
  return jvals
end

function c_lagrangian_hessian_hs71(
  userdata::Nothing,
  hvals::Vector{Float64},
  x::Vector{Float64},
  multipliers::Vector{Float64},
  objective_multiplier::Float64,
)
  # Only store the lower triangle part of the Hessian of the Lagrangian
  # Objective
  hvals[1] = objective_multiplier * (2 * x[4])                # 1,1
  hvals[2] = objective_multiplier * (x[4])                    # 2,1
  hvals[3] = 0                                                # 2,2
  hvals[4] = objective_multiplier * (x[4])                    # 3,1
  hvals[5] = 0                                                # 3,2
  hvals[6] = 0                                                # 3,3
  hvals[7] = objective_multiplier * (2 * x[1] + x[2] + x[3])  # 4,1
  hvals[8] = objective_multiplier * (x[1])                    # 4,2
  hvals[9] = objective_multiplier * (x[1])                    # 4,3
  hvals[10] = 0                                               # 4,4

  # First constraint
  hvals[2] += multipliers[1] * (x[3] * x[4])  # 2,1
  hvals[4] += multipliers[1] * (x[2] * x[4])  # 3,1
  hvals[5] += multipliers[1] * (x[1] * x[4])  # 3,2
  hvals[7] += multipliers[1] * (x[2] * x[3])  # 4,1
  hvals[8] += multipliers[1] * (x[1] * x[3])  # 4,2
  hvals[9] += multipliers[1] * (x[1] * x[2])  # 4,3

  # Second constraint
  hvals[1]  += multipliers[2] * 2  # 1,1
  hvals[3]  += multipliers[2] * 2  # 2,2
  hvals[6]  += multipliers[2] * 2  # 3,3
  hvals[10] += multipliers[2] * 2  # 4,4
  return hvals
end

function c_jacobian_operator_hs71(
  userdata::Nothing,
  Jv::Vector{Float64},
  x::Vector{Float64},
  v::Vector{Float64},
  evaluate_at_x::Bool,
)
  Jv[1] = (x[2] * x[3] * x[4]) * v[1] + (x[1] * x[3] * x[4]) * v[2] + (x[1] * x[2] * x[4]) * v[3] + (x[1] * x[2] * x[3]) * v[4]
  Jv[2] = 2 * (x[1] * v[1] + x[2] * v[2] + x[3] * v[3] + x[4] * v[4])
  return Jv
end

function c_jacobian_transposed_operator_hs71(
  userdata::Nothing,
  Jtv::Vector{Float64},
  x::Vector{Float64},
  v::Vector{Float64},
  evaluate_at_x::Bool,
)
  Jtv[1] = (x[2] * x[3] * x[4]) * v[1] + 2 * x[1] * v[2]
  Jtv[2] = (x[1] * x[3] * x[4]) * v[1] + 2 * x[2] * v[2]
  Jtv[3] = (x[1] * x[2] * x[4]) * v[1] + 2 * x[3] * v[2]
  Jtv[4] = (x[1] * x[2] * x[3]) * v[1] + 2 * x[4] * v[2]
  return Jtv
end

function c_lagrangian_hessian_operator_hs71(
  userdata::Nothing,
  Hv::Vector{Float64},
  x::Vector{Float64},
  objective_multiplier::Float64,
  multipliers::Vector{Float64},
  v::Vector{Float64},
  evaluate_at_x::Bool,
)
  # Objective
  Hv[1] = objective_multiplier * (2*x[4]*v[1] + x[4]*v[2] + x[4]*v[3] + (2*x[1] + x[2] + x[3])*v[4])
  Hv[2] = objective_multiplier * (x[4]*v[1] + 0*v[2] + 0*v[3] + x[1]*v[4])
  Hv[3] = objective_multiplier * (x[4]*v[1] + 0*v[2] + 0*v[3] + x[1]*v[4])
  Hv[4] = objective_multiplier * ((2*x[1] + x[2] + x[3])*v[1] + x[1]*v[2] + x[1]*v[3] + 0*v[4])

  # First constraint
  Hv[1] += multipliers[1] * (x[3]*x[4]*v[2] + x[2]*x[4]*v[3] + x[2]*x[3]*v[4])
  Hv[2] += multipliers[1] * (x[3]*x[4]*v[1] + x[1]*x[4]*v[3] + x[1]*x[3]*v[4])
  Hv[3] += multipliers[1] * (x[2]*x[4]*v[1] + x[1]*x[4]*v[2] + x[1]*x[2]*v[4])
  Hv[4] += multipliers[1] * (x[2]*x[3]*v[1] + x[1]*x[3]*v[2] + x[1]*x[2]*v[3])

  # Second constraint
  Hv[1] += multipliers[2] * (2*v[1])
  Hv[2] += multipliers[2] * (2*v[2])
  Hv[3] += multipliers[2] * (2*v[3])
  Hv[4] += multipliers[2] * (2*v[4])
  return Hv
end

# Sparsity pattern of the Jacobian
function sparsity_pattern_jacobian_hs71()
  nnzj = 8
  jrows = Vector{Cint}(undef, nnzj)
  jcols = Vector{Cint}(undef, nnzj)
  # Constraint (row) 1
  jrows[1] = 1
  jcols[1] = 1
  jrows[2] = 1
  jcols[2] = 2
  jrows[3] = 1
  jcols[3] = 3
  jrows[4] = 1
  jcols[4] = 4
  # Constraint (row) 2
  jrows[5] = 2
  jcols[5] = 1
  jrows[6] = 2
  jcols[6] = 2
  jrows[7] = 2
  jcols[7] = 3
  jrows[8] = 2
  jcols[8] = 4
  return jrows, jcols, nnzj
end

# Sparsity pattern of the Hessian of the Lagrangian
function sparsity_pattern_lagrangian_hessian_hs71()
  nnzh = 10
  hrows = Vector{Cint}(undef, nnzh)
  hcols = Vector{Cint}(undef, nnzh)
  # Symmetric matrix, fill the lower left triangle only
  idx = 1
  for row in 1:4
    for col in 1:row
      hrows[idx] = row
      hcols[idx] = col
      idx += 1
    end
  end
  return hrows, hcols, nnzh
end

@testset "hs71" begin
  nvar = 4
  lvar = Float64[1.0, 1.0, 1.0, 1.0]
  uvar = Float64[5.0, 5.0, 5.0, 5.0]

  ncon = 2
  lcon = Float64[25.0, 40.0]
  ucon = Float64[2.0e19, 40.0]

  jrows, jcols, nnzj = sparsity_pattern_jacobian_hs71()
  hrows, hcols, nnzh = sparsity_pattern_lagrangian_hessian_hs71()

  x0 = Float64[1.0, 5.0, 5.0, 1.0]
  y0 = zeros(Float64, ncon)

  @testset "uno_model -- uno_solver -- uno_optimize" begin
    model = uno_model(
      "NLP",
      true,
      nvar,
      ncon,
      lvar,
      uvar,
      lcon,
      ucon,
      jrows,
      jcols,
      nnzj,
      hrows,
      hcols,
      nnzh,
      c_objective_hs71,
      c_constraints_hs71,
      c_objective_gradient_hs71,
      c_jacobian_hs71,
      c_lagrangian_hessian_hs71,
      c_jacobian_operator_hs71,
      c_jacobian_transposed_operator_hs71,
      c_lagrangian_hessian_operator_hs71,
    )

    uno_set_initial_primal_iterate(model, x0)
    uno_set_initial_dual_iterate(model, y0)

    # code coverage
    uno_set_variables_lower_bounds(model, lvar)
    uno_set_variables_upper_bounds(model, uvar)
    for i = 1:nvar
      UnoSolver.uno_set_variable_lower_bound(model, i, lvar[i])
      UnoSolver.uno_set_variable_upper_bound(model, i, uvar[i])
    end
    uno_set_constraints_lower_bounds(model, lcon)
    uno_set_constraints_upper_bounds(model, ucon)
    for i = 1:ncon
      UnoSolver.uno_set_constraint_lower_bound(model, i, lcon[i])
      UnoSolver.uno_set_constraint_upper_bound(model, i, ucon[i])
    end

    solver = uno_solver()
    uno_set_solver_preset(solver, "filtersqp")
    uno_set_solver_bool_option(solver, "print_solution", true)
    uno_optimize(solver, model)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)
    hessian_model = UnoSolver.uno_get_solver_string_option(solver, "hessian_model")

    primal_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    constraint_dual_solution = Vector{Float64}(undef, ncon)
    UnoSolver.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    @test optimization_status == 0  # UNO_SUCCESS
    @test solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test primal_solution[1] ≈ 1.0000000000000000 atol = 1e-5
    @test primal_solution[2] ≈ 4.7429996418092970 atol = 1e-5
    @test primal_solution[3] ≈ 3.8211499817883077 atol = 1e-5
    @test primal_solution[4] ≈ 1.3794082897556983 atol = 1e-5
    @test solution_objective ≈ 17.014017145179164 atol = 1e-5
    @test solution_primal_feasibility ≈ 0.0 atol = 1e-5
    @test solution_stationarity ≈ 0.0 atol = 1e-5
    @test solution_complementarity ≈ 0.0 atol = 1e-5
    @test hessian_model == "exact"
  end

  @testset "L-BFGS with preset=filtersqp" begin
    model = uno_model(
      "NLP",
      true,
      nvar,
      ncon,
      lvar,
      uvar,
      lcon,
      ucon,
      jrows,
      jcols,
      nnzj,
      hrows,
      hcols,
      nnzh,
      c_objective_hs71,
      c_constraints_hs71,
      c_objective_gradient_hs71,
      c_jacobian_hs71,
      nothing,
      c_jacobian_operator_hs71,
      c_jacobian_transposed_operator_hs71,
      nothing,
    )

    uno_set_initial_primal_iterate(model, x0)
    uno_set_initial_dual_iterate(model, y0)

    solver = uno_solver()
    uno_set_solver_preset(solver, "filtersqp")
    uno_set_solver_bool_option(solver, "print_solution", true)
    uno_set_solver_string_option(solver, "hessian_model", "LBFGS")
    uno_optimize(solver, model)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)
    hessian_model = UnoSolver.uno_get_solver_string_option(solver, "hessian_model")

    primal_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    constraint_dual_solution = Vector{Float64}(undef, ncon)
    UnoSolver.uno_get_constraint_dual_solution(solver, constraint_dual_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    @test optimization_status == 0  # UNO_SUCCESS
    @test solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test primal_solution[1] ≈ 1.0000000000000000 atol = 1e-5
    @test primal_solution[2] ≈ 4.7429996418092970 atol = 1e-5
    @test primal_solution[3] ≈ 3.8211499817883077 atol = 1e-5
    @test primal_solution[4] ≈ 1.3794082897556983 atol = 1e-5
    @test solution_objective ≈ 17.014017145179164 atol = 1e-5
    @test solution_primal_feasibility ≈ 0.0 atol = 1e-5
    @test solution_stationarity ≈ 0.0 atol = 1e-5
    @test solution_complementarity ≈ 0.0 atol = 1e-5
    @test hessian_model == "LBFGS"
  end

  @testset "uno" begin
    stats = uno(
      "NLP",
      true,
      nvar,
      ncon,
      lvar,
      uvar,
      lcon,
      ucon,
      jrows,
      jcols,
      nnzj,
      hrows,
      hcols,
      nnzh,
      c_objective_hs71,
      c_constraints_hs71,
      c_objective_gradient_hs71,
      c_jacobian_hs71,
      c_lagrangian_hessian_hs71,
      c_jacobian_operator_hs71,
      c_jacobian_transposed_operator_hs71,
      c_lagrangian_hessian_operator_hs71,
      nothing,
      'L',
      1,
      x0,
      y0;
      preset="filtersqp",
      print_solution=true,
    )

    @test stats.optimization_status == 0  # UNO_SUCCESS
    @test stats.solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test stats.solution[1] ≈ 1.0000000000000000 atol = 1e-5
    @test stats.solution[2] ≈ 4.7429996418092970 atol = 1e-5
    @test stats.solution[3] ≈ 3.8211499817883077 atol = 1e-5
    @test stats.solution[4] ≈ 1.3794082897556983 atol = 1e-5
    @test stats.objective ≈ 17.014017145179164 atol = 1e-5
    @test stats.primal_feas ≈ 0.0 atol = 1e-5
    @test stats.dual_feas ≈ 0.0 atol = 1e-5
    @test stats.complementarity_feas ≈ 0.0 atol = 1e-5
  end
end

# hs5
# min sin(x1 + x2) + (x1 - x2)^2 - 1.5 * x1 + 2.5 * x2 + 1
# st  -1.5 <= x1 <= -3, 4 <= x2 <= 3
# Start at (0, 0)
# End at (-0.547..., -1.547...)

function c_objective_hs5(userdata::Nothing, x::Vector{Float64})
  f = sin(x[1] + x[2]) + (x[1] - x[2])^2 - 1.5 * x[1] + 2.5 * x[2] + 1.0
  return f
end

function c_objective_gradient_hs5(userdata::Nothing, g::Vector{Float64}, x::Vector{Float64})
  g[1] = cos(x[1] + x[2]) + 2 * (x[1] - x[2]) - 1.5
  g[2] = cos(x[1] + x[2]) - 2 * (x[1] - x[2]) + 2.5
  return g
end

function c_objective_hessian_hs5(
  userdata::Nothing,
  hvals::Vector{Float64},
  x::Vector{Float64},
  multipliers::Vector{Float64}, # ignored
  objective_multiplier::Float64,
)
  # Only store the upper triangle part of the Hessian of the objective
  # Objective
  hvals[1] = objective_multiplier * (-sin(x[1] + x[2]) + 2)  # 1,1
  hvals[2] = objective_multiplier * (-sin(x[1] + x[2]) - 2)  # 1,2
  hvals[3] = objective_multiplier * (-sin(x[1] + x[2]) + 2)  # 2,2
  return hvals
end

function c_objective_hessian_operator_hs5(
  userdata::Nothing,
  Hv::Vector{Float64},
  x::Vector{Float64},
  objective_multiplier::Float64,
  multipliers::Vector{Float64}, # ignored
  v::Vector{Float64},
  evaluate_at_x::Bool,
)
  # Objective
  Hv[1] = (-sin(x[1] + x[2]) * (v[1] + v[2]) + 2 * (v[1] - v[2])) * objective_multiplier
  Hv[2] = (-sin(x[1] + x[2]) * (v[1] + v[2]) + 2 * (v[2] - v[1])) * objective_multiplier
  return Hv
end

# Sparsity pattern of the Hessian of the objective
function sparsity_pattern_objective_hessian_hs5()
  nnzh = 3
  hrows = Vector{Cint}(undef, nnzh)
  hcols = Vector{Cint}(undef, nnzh)
  # Symmetric matrix, fill the upper right triangle only
  idx = 1
  for row in 1:2
    for col in row:2
      hrows[idx] = row
      hcols[idx] = col
      idx += 1
    end
  end
  return hrows, hcols, nnzh
end

@testset "hs5" begin
  nvar = 2
  lvar = Float64[-1.5, -3.0]
  uvar = Float64[4.0, 3.0]

  hrows, hcols, nnzh = sparsity_pattern_objective_hessian_hs5()

  x0 = Float64[0.0, 0.0]

  @testset "uno_model -- uno_solver -- uno_optimize" begin
    model = uno_model(
      "NLP",
      true,
      nvar,
      lvar,
      uvar,
      hrows,
      hcols,
      nnzh,
      c_objective_hs5,
      c_objective_gradient_hs5,
      c_objective_hessian_hs5,
      c_objective_hessian_operator_hs5,
      nothing,
      'U',
      1,
    )

    uno_set_initial_primal_iterate(model, x0)

    # code coverage
    uno_set_variables_lower_bounds(model, lvar)
    uno_set_variables_upper_bounds(model, uvar)
    for i = 1:nvar
      UnoSolver.uno_set_variable_lower_bound(model, i, lvar[i])
      UnoSolver.uno_set_variable_upper_bound(model, i, uvar[i])
    end

    solver = uno_solver()
    uno_set_solver_preset(solver, "filtersqp")
    uno_set_solver_bool_option(solver, "print_solution", true)
    uno_optimize(solver, model)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)
    hessian_model = UnoSolver.uno_get_solver_string_option(solver, "hessian_model")

    primal_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    @test optimization_status == 0  # UNO_SUCCESS
    @test solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test primal_solution[1] ≈ -0.5471975484809134 atol = 1e-5
    @test primal_solution[2] ≈ -1.5471975484809133 atol = 1e-5
    @test solution_objective ≈ -1.9132229549810362 atol = 1e-5
    @test solution_primal_feasibility ≈ 0.0 atol = 1e-5
    @test solution_stationarity ≈ 0.0 atol = 1e-5
    @test solution_complementarity ≈ 0.0 atol = 1e-5
    @test hessian_model == "exact"
  end

  @testset "L-BFGS with preset=filtersqp" begin
    model = uno_model(
      "NLP",
      true,
      nvar,
      lvar,
      uvar,
      hrows,
      hcols,
      nnzh,
      c_objective_hs5,
      c_objective_gradient_hs5,
      nothing,
      nothing,
      nothing,
      'U',
      1,
    )

    uno_set_initial_primal_iterate(model, x0)

    solver = uno_solver()
    uno_set_solver_preset(solver, "filtersqp")
    uno_set_solver_bool_option(solver, "print_solution", true)
    uno_set_solver_string_option(solver, "hessian_model", "LBFGS")
    uno_optimize(solver, model)

    optimization_status = UnoSolver.uno_get_optimization_status(solver)
    solution_status = UnoSolver.uno_get_solution_status(solver)
    solution_objective = UnoSolver.uno_get_solution_objective(solver)
    solution_primal_feasibility = UnoSolver.uno_get_solution_primal_feasibility(solver)
    solution_stationarity = UnoSolver.uno_get_solution_stationarity(solver)
    solution_complementarity = UnoSolver.uno_get_solution_complementarity(solver)
    hessian_model = UnoSolver.uno_get_solver_string_option(solver, "hessian_model")

    primal_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_primal_solution(solver, primal_solution)

    lower_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_lower_bound_dual_solution(solver, lower_bound_dual_solution)

    upper_bound_dual_solution = Vector{Float64}(undef, nvar)
    UnoSolver.uno_get_upper_bound_dual_solution(solver, upper_bound_dual_solution)

    @test optimization_status == 0  # UNO_SUCCESS
    @test solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test primal_solution[1] ≈ -0.5471975484809134 atol = 1e-5
    @test primal_solution[2] ≈ -1.5471975484809133 atol = 1e-5
    @test solution_objective ≈ -1.9132229549810362 atol = 1e-5
    @test solution_primal_feasibility ≈ 0.0 atol = 1e-5
    @test solution_stationarity ≈ 0.0 atol = 1e-5
    @test solution_complementarity ≈ 0.0 atol = 1e-5
    @test hessian_model == "LBFGS"
  end

  @testset "uno" begin
    stats = uno(
      "NLP",
      true,
      nvar,
      lvar,
      uvar,
      hrows,
      hcols,
      nnzh,
      c_objective_hs5,
      c_objective_gradient_hs5,
      c_objective_hessian_hs5,
      c_objective_hessian_operator_hs5,
      nothing,
      'U',
      1,
      x0,
      preset="filtersqp",
      print_solution=true,
    )

    @test stats.optimization_status == 0  # UNO_SUCCESS
    @test stats.solution_status == 1      # UNO_FEASIBLE_KKT_POINT
    @test stats.solution[1] ≈ -0.5471975484809134 atol = 1e-5
    @test stats.solution[2] ≈ -1.5471975484809133 atol = 1e-5
    @test stats.objective ≈ -1.9132229549810362 atol = 1e-5
    @test stats.primal_feas ≈ 0.0 atol = 1e-5
    @test stats.dual_feas ≈ 0.0 atol = 1e-5
    @test stats.complementarity_feas ≈ 0.0 atol = 1e-5
  end
end
