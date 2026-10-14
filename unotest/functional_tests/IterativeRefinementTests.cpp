// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Iterative refinement on the uncondensed Newton system of the primal-dual interior-point method, and singularity
// handling in the inertia correction strategy

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/hessian_models/ExactHessian.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/InteriorPointParameters.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "ingredients/inertia_correction_strategies/PrimalDualInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/EQPSolver.hpp"
#include "ingredients/subproblem_solvers/LinearSystem.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "symbolic/IntegerRange.hpp"
#include "tools/Infinity.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   // min  (x0 - 1)^2 + (x1 - 2)^2 + x0 x1
   // s.t. x0 + x1^2 = 1
   //      x0 - x1 <= 3   (slacked)
   //      0 <= x0 <= 10, x1 >= -5
   class BoundedModel: public Model {
   public:
      BoundedModel(): Model("bounded", 2, 2, 1., -1., 0) { }

      ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      bool has_jacobian_operator() const override { return false; }
      bool has_jacobian_transposed_operator() const override { return false; }
      bool has_hessian_operator() const override { return false; }
      bool has_hessian_matrix() const override { return true; }

      double evaluate_objective(const Vector<double>& x) const override {
         return std::pow(x[0] - 1., 2) + std::pow(x[1] - 2., 2) + x[0] * x[1];
      }
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         constraints[0] = x[0] + x[1] * x[1];
         constraints[1] = x[0] - x[1];
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         gradient[0] = 2. * (x[0] - 1.) + x[1];
         gradient[1] = 2. * (x[1] - 2.) + x[0];
      }
      View<const uno_int> get_jacobian_row_indices() const override { return view(this->jacobian_rows.data(), 4); }
      View<const uno_int> get_jacobian_column_indices() const override { return view(this->jacobian_columns.data(), 4); }
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         jacobian_values[0] = 1.;
         jacobian_values[1] = 2. * x[1];
         jacobian_values[2] = 1.;
         jacobian_values[3] = -1.;
      }
      // lower triangle: (0, 0), (1, 1), (1, 0)
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         const uno_int rows[] = {0, 1, 1};
         const uno_int columns[] = {0, 1, 0};
         for (size_t index = 0; index < 3; ++index) {
            row_indices[index] = rows[index] + solver_indexing;
            column_indices[index] = columns[index] + solver_indexing;
         }
      }
      // Lagrangian rho f(x) - y^T c(x)
      void evaluate_lagrangian_hessian(const Vector<double>& /*x*/, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         hessian_values[0] = 2. * objective_multiplier;
         hessian_values[1] = 2. * objective_multiplier - 2. * multipliers[0];
         hessian_values[2] = objective_multiplier;
      }
      void compute_jacobian_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_jacobian_transposed_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_hessian_vector_product(View<const double>, View<const double>, double, const Vector<double>&,
            View<double>) const override {
         throw std::runtime_error("not implemented");
      }

      const std::vector<double>& get_variables_lower_bounds() const override { return this->variables_lower_bounds; }
      const std::vector<double>& get_variables_upper_bounds() const override { return this->variables_upper_bounds; }
      const std::vector<double>& get_constraints_lower_bounds() const override { return this->constraints_lower_bounds; }
      const std::vector<double>& get_constraints_upper_bounds() const override { return this->constraints_upper_bounds; }
      const std::vector<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
      const std::vector<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
      const Collection<size_t>& get_linear_constraints() const override { return this->linear_constraints; }
      const Collection<size_t>& get_nonlinear_constraints() const override { return this->nonlinear_constraints; }

      void initial_primal_point(Vector<double>& x) const override {
         x[0] = 0.5;
         x[1] = 0.3;
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      size_t number_jacobian_nonzeros() const override { return 4; }
      size_t number_hessian_nonzeros() const override { return 3; }

      size_t number_model_objective_evaluations() const override { return 0; }
      size_t number_model_constraints_evaluations() const override { return 0; }
      size_t number_model_objective_gradient_evaluations() const override { return 0; }
      size_t number_model_jacobian_evaluations() const override { return 0; }
      size_t number_model_hessian_evaluations() const override { return 0; }
      void reset_number_evaluations() const override { }

   private:
      const std::vector<uno_int> jacobian_rows{0, 0, 1, 1};
      const std::vector<uno_int> jacobian_columns{0, 1, 0, 1};
      const std::vector<double> variables_lower_bounds{0., -5.};
      const std::vector<double> variables_upper_bounds{10., INFINITY};
      const std::vector<double> constraints_lower_bounds{1., -INFINITY};
      const std::vector<double> constraints_upper_bounds{1., 3.};
      const std::vector<size_t> equality_constraints{0};
      const std::vector<size_t> inequality_constraints{1};
      const IntegerRange linear_constraints{1, 2};
      const IntegerRange nonlinear_constraints{0, 1};
   };

   Options make_options(const std::string& linear_solver, size_t min_steps, size_t max_steps) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_string("linear_solver", linear_solver);
      options.set_integer("iterative_refinement_min_steps", static_cast<uno_int>(min_steps));
      options.set_integer("iterative_refinement_max_steps", static_cast<uno_int>(max_steps));
      return options;
   }

   // options with which iterative refinement always fails: the solver's failure handling is triggered
   Options make_failing_options(const std::string& linear_solver) {
      Options options = make_options(linear_solver, 1, 10);
      options.set_double("iterative_refinement_residual_ratio_max", 0.);
      options.set_double("iterative_refinement_residual_ratio_singular", 0.);
      return options;
   }

   bool supports_pivoting_tolerance_increase(const std::string& linear_solver) {
      return linear_solver == "MA27" || linear_solver == "MA57" || linear_solver == "MUMPS";
   }

   // primal-dual inertia correction that counts the inertia-checked refactorizations
   class CountingInertiaCorrection: public PrimalDualInertiaCorrection {
   public:
      using PrimalDualInertiaCorrection::PrimalDualInertiaCorrection;
      size_t number_refactorizations{0};

      void refactorize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem, double dual_regularization_parameter,
            const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
            BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) override {
         ++this->number_refactorizations;
         PrimalDualInertiaCorrection::refactorize_augmented_matrix(statistics, subproblem, dual_regularization_parameter,
            expected_inertia, linear_solver, block_augmented_matrix, primal_regularization_variables);
      }
   };

   // the primal-dual IPM ingredients at the initial iterate, with a strictly positive barrier parameter
   struct BarrierSetting {
      BoundedModel model{};
      Options options;
      OptimizationProblem problem{model};
      InteriorPointParameters parameters{0.99, 1e10, 0.25, 10., 1e-2, 1e-2, 1e-5, 1., 1e-12};
      double barrier_parameter{0.1};
      PrimalDualInteriorPointProblem barrier_problem{problem, parameters, barrier_parameter};
      ExactHessian hessian_model{model};
      CountingInertiaCorrection inertia_correction;
      Subproblem subproblem{barrier_problem, hessian_model, inertia_correction};
      Statistics statistics{};
      Evaluations evaluations{model};
      Iterate iterate{};

      explicit BarrierSetting(const Options& options): options(options), inertia_correction(this->options) {
         this->hessian_model.initialize_statistics(this->statistics);
         this->inertia_correction.initialize_statistics(this->statistics);
         this->iterate = Iterate(this->model.number_variables, this->model.number_constraints);
         this->model.initial_primal_point(this->iterate.primals);
         this->barrier_problem.create_initial_iterate(this->iterate, this->evaluations);
         this->iterate.multipliers.constraints.fill(0.5);
      }

      const Direction& solve(EQPSolver<NoLowRankCorrection>& solver) {
         solver.initialize_memory(this->subproblem);
         const WarmstartInformation warmstart_information{};
         const Vector<double> initial_point(this->barrier_problem.number_variables);
         return solver.solve(this->statistics, this->subproblem, this->iterate, Inf, initial_point,
            this->evaluations, warmstart_information);
      }
   };

   // iterative refinement is tested with every available direct linear solver
   class IterativeRefinement: public testing::TestWithParam<std::string> { };

   // access to the linear solver and the refinement tolerance (to simulate failing, then accurate iterations)
   class InspectableEQPSolver: public EQPSolver<NoLowRankCorrection> {
   public:
      using EQPSolver<NoLowRankCorrection>::EQPSolver;
      DirectSymmetricIndefiniteLinearSolver<double>* linear_solver_pointer() { return this->linear_solver.get(); }
      void set_residual_ratio_max(double value) {
         this->iterative_refinement.residual_ratio_max = value;
         this->iterative_refinement.residual_ratio_singular = std::max(value, this->iterative_refinement.residual_ratio_singular);
      }
   };
} // namespace

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(IterativeRefinement);
INSTANTIATE_TEST_SUITE_P(LinearSolvers, IterativeRefinement,
   testing::ValuesIn(SymmetricIndefiniteLinearSolverFactory::available_solvers()),
   [](const testing::TestParamInfo<std::string>& info) { return info.param; });

// condensing the uncondensed Newton RHS gives back the augmented RHS
TEST_P(IterativeRefinement, CondensedNewtonRHS) {
   BarrierSetting setting(make_options(GetParam(), 1, 10));
   const size_t augmented_dimension = setting.barrier_problem.number_variables + setting.barrier_problem.number_constraints;
   Vector<double> augmented_rhs(augmented_dimension);
   setting.subproblem.assemble_augmented_rhs(setting.iterate, setting.evaluations, augmented_rhs);
   Vector<double> uncondensed_rhs(setting.barrier_problem.uncondensed_system_dimension());
   setting.barrier_problem.uncondense_rhs(setting.iterate, augmented_rhs, uncondensed_rhs);
   Vector<double> condensed_rhs(augmented_dimension);
   setting.barrier_problem.condense_rhs(setting.iterate, uncondensed_rhs, condensed_rhs);
   for (size_t index: Range(augmented_dimension)) {
      EXPECT_NEAR(condensed_rhs[index], augmented_rhs[index], 1e-12 * (1. + std::abs(augmented_rhs[index])));
   }
}

// the refined direction solves the uncondensed Newton system. The residuals are computed independently of the
// iterative refinement machinery, from the derivatives of the model:
//  stationarity:    W dx - J^T dy - dz_L - dz_U = -(∇f + damping - J^T y - z_L - z_U)
//  constraints:     J dx = -c(x)
//  complementarity: z_L dx + (x - l) dz_L = mu - (x - l) z_L   (idem for the upper bounds)
TEST_P(IterativeRefinement, DirectionSolvesUncondensedSystem) {
   BarrierSetting setting(make_options(GetParam(), 1, 10));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   const Direction& direction = setting.solve(solver);
   EXPECT_LE(solver.get_residual_ratio(), 1e-10);
   // the matrix is not regularized at this iterate (the residuals below assume it)
   ASSERT_EQ(setting.inertia_correction.get_primal_regularization_factor(), 0.);

   // barrier problem: variables (x0, x1, s), constraints c0 = x0 + x1^2 - 1 = 0 and c1 = x0 - x1 - s = 0
   const auto& x = setting.iterate.primals;
   const auto& y = setting.iterate.multipliers.constraints;
   const auto& z_L = setting.iterate.multipliers.lower_bounds;
   const auto& z_U = setting.iterate.multipliers.upper_bounds;
   const auto& dx = direction.primals;
   const auto& dy = direction.multipliers.constraints;
   const auto& dz_L = direction.multipliers.lower_bounds;
   const auto& dz_U = direction.multipliers.upper_bounds;
   const double mu = setting.barrier_parameter;
   const double damping = setting.parameters.damping_factor * mu;
   ASSERT_EQ(x.size(), 3);

   const double lower_bounds[3] = {0., -5., -INFINITY};
   const double upper_bounds[3] = {10., INFINITY, 3.};
   const double W[3][3] = {{2., 1., 0.}, {1., 2. - 2. * y[0], 0.}, {0., 0., 0.}};
   const double J[2][3] = {{1., 2. * x[1], 0.}, {1., -1., -1.}};
   const double gradient[3] = {2. * (x[0] - 1.) + x[1], 2. * (x[1] - 2.) + x[0], 0.};
   // damping of the variables bounded on one side only: x1 (lower), s (upper)
   const double damping_terms[3] = {0., damping, -damping};
   const double constraints[2] = {x[0] + x[1] * x[1] - 1., x[0] - x[1] - x[2]};

   const auto expect_near = [](double lhs, double rhs) {
      EXPECT_NEAR(lhs, rhs, 1e-9 * (1. + std::abs(rhs)));
   };
   for (size_t i: Range(3)) {
      double lhs = -dz_L[i] - dz_U[i];
      double rhs = -(gradient[i] + damping_terms[i] - z_L[i] - z_U[i]);
      for (size_t j: Range(3)) {
         lhs += W[i][j] * dx[j];
      }
      for (size_t k: Range(2)) {
         lhs -= J[k][i] * dy[k];
         rhs += J[k][i] * y[k];
      }
      expect_near(lhs, rhs);
   }
   for (size_t k: Range(2)) {
      double lhs = 0.;
      for (size_t j: Range(3)) {
         lhs += J[k][j] * dx[j];
      }
      expect_near(lhs, -constraints[k]);
   }
   for (size_t i: Range(3)) {
      if (std::isfinite(lower_bounds[i])) {
         const double distance = x[i] - lower_bounds[i];
         expect_near(z_L[i] * dx[i] + distance * dz_L[i], mu - distance * z_L[i]);
      }
      if (std::isfinite(upper_bounds[i])) {
         const double distance = x[i] - upper_bounds[i];
         expect_near(z_U[i] * dx[i] + distance * dz_U[i], mu - distance * z_U[i]);
      }
   }
}

// default: no refinement step when the solution is already accurate
TEST_P(IterativeRefinement, DefaultSkipsRefinementWhenAccurate) {
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   options.set_string("linear_solver", GetParam());
   BarrierSetting setting(options);
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 0);
   EXPECT_LE(solver.get_residual_ratio(), 1e-10);
}

// IPOPT's default (min_steps = 1): one refinement step per solve
TEST_P(IterativeRefinement, IpoptMinimumPerformsOneStep) {
   BarrierSetting setting(make_options(GetParam(), 1, 10));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 1);
}

// max_steps = 0 disables iterative refinement altogether (min_steps is then irrelevant)
TEST_P(IterativeRefinement, DisabledWhenMaxIsZero) {
   BarrierSetting setting(make_options(GetParam(), 1, 0));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   const Direction& direction = setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 0);
   EXPECT_TRUE(std::isnan(solver.get_residual_ratio()));
   EXPECT_FALSE(solver.is_matrix_deemed_singular());
   EXPECT_TRUE(std::isfinite(direction.norm));
}

TEST_P(IterativeRefinement, MinimumNumberOfSteps) {
   BarrierSetting setting(make_options(GetParam(), 3, 3));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 3);
}

TEST_P(IterativeRefinement, InconsistentOptions) {
   EXPECT_THROW(EQPSolver<NoLowRankCorrection> solver(make_options(GetParam(), 2, 1)), std::invalid_argument);
   Options options = make_options(GetParam(), 1, 10);
   options.set_double("iterative_refinement_improvement_factor", 0.);
   EXPECT_THROW(EQPSolver<NoLowRankCorrection> solver(options), std::invalid_argument);
   options.set_double("iterative_refinement_improvement_factor", 1.5);
   EXPECT_THROW(EQPSolver<NoLowRankCorrection> solver(options), std::invalid_argument);
   options = make_options(GetParam(), 1, 10);
   options.set_double("iterative_refinement_residual_ratio_max", -1.);
   EXPECT_THROW(EQPSolver<NoLowRankCorrection> solver(options), std::invalid_argument);
}

// failure path inside the solver: the pivoting tolerance is increased (if supported) and the matrix refactorized with an
// inertia check, then the matrix is deemed singular and regularized. The direction is still computed
TEST_P(IterativeRefinement, FailureIncreasesPivotingToleranceThenRegularizes) {
   BarrierSetting setting(make_failing_options(GetParam()));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   const Direction& direction = setting.solve(solver);
   ASSERT_GT(solver.get_residual_ratio(), 0.) << "the residual ratio is exactly zero: the failure path cannot be tested";
   const size_t expected_increases = supports_pivoting_tolerance_increase(GetParam()) ? 1 : 0;
   // the refactorization with the increased pivoting tolerance goes through the inertia correction strategy
   EXPECT_EQ(setting.inertia_correction.number_refactorizations, expected_increases);
   EXPECT_TRUE(solver.is_matrix_deemed_singular());
   // pretend singular: the dual regularization is added first. The dual regularization block ends the matrix values
   const auto& linear_system = static_cast<const LinearSystem&>(solver.get_workspace());
   EXPECT_LT(linear_system.matrix_values[linear_system.number_nonzeros - 1], 0.);
   EXPECT_NE(direction.status, SubproblemStatus::INFEASIBLE);
   EXPECT_TRUE(std::isfinite(direction.norm));
}

// second-order corrections never modify the factorization, even when refinement fails
TEST_P(IterativeRefinement, SecondOrderCorrectionDoesNotModifyFactorization) {
   BarrierSetting setting(make_failing_options(GetParam()));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   const double primal_regularization = setting.inertia_correction.get_primal_regularization_factor();

   Iterate trial_iterate = setting.iterate;
   Evaluations trial_evaluations{setting.model};
   solver.initialize_second_order_corrections(setting.subproblem, setting.iterate, trial_iterate, setting.evaluations,
      trial_evaluations);
   const Direction& direction = solver.compute_second_order_correction(setting.subproblem, setting.iterate);
   EXPECT_FALSE(solver.is_matrix_deemed_singular());
   EXPECT_EQ(setting.inertia_correction.get_primal_regularization_factor(), primal_regularization);
   for (size_t variable_index: Range(setting.barrier_problem.number_variables)) {
      EXPECT_TRUE(std::isfinite(direction.primals[variable_index]));
   }
}

// a matrix deemed singular is first given a dual regularization, then an increasing primal regularization
TEST_P(IterativeRefinement, SingularMatrixIsRegularized) {
   BarrierSetting setting(make_options(GetParam(), 1, 10));
   const auto linear_solver = SymmetricIndefiniteLinearSolverFactory::create(setting.options);
   auto& linear_system = linear_solver->get_linear_system();
   linear_system.initialize_augmented_system(setting.subproblem);
   linear_solver->initialize_memory();
   auto block_augmented_matrix = setting.subproblem.compute_block_augmented_matrix(linear_system);
   setting.subproblem.evaluate_lagrangian_hessian(setting.statistics, setting.iterate, block_augmented_matrix.hessian);
   setting.subproblem.evaluate_jacobian(setting.iterate, block_augmented_matrix.jacobian, setting.evaluations);
   linear_solver->do_symbolic_analysis();
   const double dual_regularization_parameter = setting.subproblem.dual_regularization_factor();
   setting.subproblem.regularize_augmented_matrix(setting.statistics, dual_regularization_parameter, *linear_solver,
      block_augmented_matrix);
   // the inertia is correct without regularization
   EXPECT_EQ(setting.inertia_correction.get_primal_regularization_factor(), 0.);
   EXPECT_EQ(block_augmented_matrix.dual_inertia_correction[0], 0.);

   // 1st singularity: dual regularization
   ASSERT_TRUE(setting.subproblem.regularize_singular_augmented_matrix(setting.statistics, dual_regularization_parameter,
      *linear_solver, block_augmented_matrix));
   const double expected_dual_regularization = setting.options.get_double("dual_regularization_fraction") * dual_regularization_parameter;
   EXPECT_EQ(setting.inertia_correction.get_primal_regularization_factor(), 0.);
   EXPECT_DOUBLE_EQ(block_augmented_matrix.dual_inertia_correction[0], -expected_dual_regularization);
   EXPECT_FALSE(linear_solver->matrix_is_singular());

   // 2nd singularity: primal regularization
   ASSERT_TRUE(setting.subproblem.regularize_singular_augmented_matrix(setting.statistics, dual_regularization_parameter,
      *linear_solver, block_augmented_matrix));
   EXPECT_DOUBLE_EQ(setting.inertia_correction.get_primal_regularization_factor(),
      setting.options.get_double("primal_regularization_initial_factor"));
   EXPECT_DOUBLE_EQ(block_augmented_matrix.primal_inertia_correction[0], setting.options.get_double("primal_regularization_initial_factor"));

   // refactorization (e.g. after a pivoting tolerance increase) keeps the current regularization if the inertia is correct
   setting.subproblem.refactorize_augmented_matrix(setting.statistics, dual_regularization_parameter, *linear_solver,
      block_augmented_matrix);
   EXPECT_DOUBLE_EQ(setting.inertia_correction.get_primal_regularization_factor(),
      setting.options.get_double("primal_regularization_initial_factor"));
   EXPECT_DOUBLE_EQ(block_augmented_matrix.dual_inertia_correction[0], -expected_dual_regularization);
   EXPECT_FALSE(linear_solver->matrix_is_singular());
}
