// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Iterative refinement on the uncondensed Newton system of the primal-dual interior-point method, and singularity
// handling in the inertia correction strategy

#ifdef HAS_MUMPS

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
#include "ingredients/subproblem_solvers/MUMPS/MUMPSSolver.hpp"
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

   Options make_options(size_t min_steps, size_t max_steps) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_string("linear_solver", "MUMPS");
      options.set_integer("iterative_refinement_min_steps", static_cast<uno_int>(min_steps));
      options.set_integer("iterative_refinement_max_steps", static_cast<uno_int>(max_steps));
      return options;
   }

   // the primal-dual IPM ingredients at the initial iterate, with a strictly positive barrier parameter
   struct BarrierSetting {
      BoundedModel model{};
      Options options;
      OptimizationProblem problem{model};
      InteriorPointParameters parameters{0.99, 1e10, 0.25, 10., 1e-2, 1e-2, 1e-5, 1., 1e-12};
      double barrier_parameter{0.1};
      PrimalDualInteriorPointProblem barrier_problem{problem, parameters, barrier_parameter};
      ExactHessian hessian_model{model};
      PrimalDualInertiaCorrection inertia_correction;
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
} // namespace

// condensing the uncondensed Newton RHS gives back the augmented RHS
TEST(IterativeRefinement, CondensedNewtonRHS) {
   BarrierSetting setting(make_options(1, 10));
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

// the refined direction solves the uncondensed system, incl. the complementarity rows
TEST(IterativeRefinement, DirectionSolvesUncondensedSystem) {
   BarrierSetting setting(make_options(1, 10));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   const Direction& direction = setting.solve(solver);
   EXPECT_LE(solver.get_residual_ratio(), 1e-10);
   // Z_L dx + S_L dz_L = mu - S_L z_L
   const auto& lower_bounds = setting.model.get_variables_lower_bounds();
   for (size_t variable_index: Range(setting.model.number_variables)) {
      const double distance = setting.iterate.primals[variable_index] - lower_bounds[variable_index];
      const double lhs = setting.iterate.multipliers.lower_bounds[variable_index] * direction.primals[variable_index] +
         distance * direction.multipliers.lower_bounds[variable_index];
      EXPECT_NEAR(lhs, setting.barrier_parameter - distance * setting.iterate.multipliers.lower_bounds[variable_index], 1e-12);
   }
}

// default (IPOPT): at least one refinement step per solve
TEST(IterativeRefinement, DefaultPerformsOneStep) {
   BarrierSetting setting(make_options(1, 10));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 1);
}

// unlike IPOPT, max_steps = 0 means no refinement step at all
TEST(IterativeRefinement, NoStepWhenMaxIsZero) {
   BarrierSetting setting(make_options(0, 0));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   const Direction& direction = setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 0);
   EXPECT_LE(solver.get_residual_ratio(), 1e-10);
   EXPECT_TRUE(std::isfinite(direction.norm));
}

TEST(IterativeRefinement, MinimumNumberOfSteps) {
   BarrierSetting setting(make_options(3, 3));
   EQPSolver<NoLowRankCorrection> solver(setting.options);
   (void)setting.solve(solver);
   EXPECT_EQ(solver.get_number_refinement_steps(), 3);
}

TEST(IterativeRefinement, InconsistentOptions) {
   const Options options = make_options(2, 1);
   EXPECT_THROW(EQPSolver<NoLowRankCorrection> solver(options), std::invalid_argument);
}

// a matrix deemed singular is first given a dual regularization, then an increasing primal regularization
TEST(IterativeRefinement, SingularMatrixIsRegularized) {
   BarrierSetting setting(make_options(1, 10));
   MUMPSSolver linear_solver;
   auto& linear_system = linear_solver.get_linear_system();
   linear_system.initialize_augmented_system(setting.subproblem);
   linear_solver.initialize_memory();
   auto block_augmented_matrix = setting.subproblem.compute_block_augmented_matrix(linear_system);
   setting.subproblem.evaluate_lagrangian_hessian(setting.statistics, setting.iterate, block_augmented_matrix.hessian);
   setting.subproblem.evaluate_jacobian(setting.iterate, block_augmented_matrix.jacobian, setting.evaluations);
   linear_solver.do_symbolic_analysis();
   const double dual_regularization_parameter = setting.subproblem.dual_regularization_factor();
   setting.subproblem.regularize_augmented_matrix(setting.statistics, dual_regularization_parameter, linear_solver,
      block_augmented_matrix);
   // the inertia is correct without regularization
   EXPECT_EQ(setting.inertia_correction.get_primal_regularization_factor(), 0.);
   EXPECT_EQ(block_augmented_matrix.dual_inertia_correction[0], 0.);

   // 1st singularity: dual regularization
   ASSERT_TRUE(setting.subproblem.regularize_singular_augmented_matrix(setting.statistics, dual_regularization_parameter,
      linear_solver, block_augmented_matrix));
   const double expected_dual_regularization = setting.options.get_double("dual_regularization_fraction") * dual_regularization_parameter;
   EXPECT_EQ(setting.inertia_correction.get_primal_regularization_factor(), 0.);
   EXPECT_DOUBLE_EQ(block_augmented_matrix.dual_inertia_correction[0], -expected_dual_regularization);
   EXPECT_FALSE(linear_solver.matrix_is_singular());

   // 2nd singularity: primal regularization
   ASSERT_TRUE(setting.subproblem.regularize_singular_augmented_matrix(setting.statistics, dual_regularization_parameter,
      linear_solver, block_augmented_matrix));
   EXPECT_DOUBLE_EQ(setting.inertia_correction.get_primal_regularization_factor(),
      setting.options.get_double("primal_regularization_initial_factor"));
   EXPECT_DOUBLE_EQ(block_augmented_matrix.primal_inertia_correction[0], setting.options.get_double("primal_regularization_initial_factor"));
}

#endif // HAS_MUMPS
