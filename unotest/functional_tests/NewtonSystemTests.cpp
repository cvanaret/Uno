// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Symbolic Newton system: the matrices assembled from the symbols coincide (pattern and values, bitwise) with the
// augmented matrices of the subproblem (unconstrained/equality-constrained problem and primal-dual interior-point problem),
// and the operator products of the symbols coincide with the products with their assembled entries

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/hessian_models/ExactHessian.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/InteriorPointParameters.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "ingredients/inertia_correction_strategies/PrimalDualInertiaCorrection.hpp"
#include "ingredients/inertia_correction_strategies/PrimalInertiaCorrection.hpp"
#include "ingredients/subproblem/NewtonSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/COOLinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/Indexing.hpp"
#include "model/Model.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "symbolic/IntegerRange.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   // min  (x0 - 1)^2 + (x1 - 2)^2 + x0 x1 + x2^2 + (x3 - 1)^4 + x4^2
   // s.t. x0 + x1^2 + x2 = 1            (equality)
   //      x0 - x1 <= 3                   (inequality: slack with an upper bound)
   //      -1 <= x2 + x3 <= 2             (range: slack with two bounds)
   //      0 <= x0 <= 10, x1 >= -5, x2 <= 4, x3 free, x4 = 0.5 (fixed)
   class BarrierTestModel: public Model {
   public:
      BarrierTestModel(): Model("barrier_test", 5, 3, 1., -1., 0) { }

      ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      bool has_jacobian_operator() const override { return false; }
      bool has_jacobian_transposed_operator() const override { return false; }
      bool has_hessian_operator() const override { return false; }
      bool has_hessian_matrix() const override { return true; }

      double evaluate_objective(const Vector<double>& x) const override {
         return std::pow(x[0] - 1., 2) + std::pow(x[1] - 2., 2) + x[0] * x[1] + x[2] * x[2] + std::pow(x[3] - 1., 4) + x[4] * x[4];
      }
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         constraints[0] = x[0] + x[1] * x[1] + x[2];
         constraints[1] = x[0] - x[1];
         constraints[2] = x[2] + x[3];
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         gradient[0] = 2. * (x[0] - 1.) + x[1];
         gradient[1] = 2. * (x[1] - 2.) + x[0];
         gradient[2] = 2. * x[2];
         gradient[3] = 4. * std::pow(x[3] - 1., 3);
         gradient[4] = 2. * x[4];
      }
      View<const uno_int> get_jacobian_row_indices() const override { return view(this->jacobian_rows.data(), 7); }
      View<const uno_int> get_jacobian_column_indices() const override { return view(this->jacobian_columns.data(), 7); }
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         jacobian_values[0] = 1.;
         jacobian_values[1] = 2. * x[1];
         jacobian_values[2] = 1.;
         jacobian_values[3] = 1.;
         jacobian_values[4] = -1.;
         jacobian_values[5] = 1.;
         jacobian_values[6] = 1.;
      }
      // lower triangle: (0, 0), (1, 1), (1, 0), (2, 2), (3, 3), (4, 4)
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         const uno_int rows[] = {0, 1, 1, 2, 3, 4};
         const uno_int columns[] = {0, 1, 0, 2, 3, 4};
         for (size_t index = 0; index < 6; ++index) {
            row_indices[index] = rows[index] + solver_indexing;
            column_indices[index] = columns[index] + solver_indexing;
         }
      }
      // Lagrangian rho f(x) - y^T c(x)
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         hessian_values[0] = 2. * objective_multiplier;
         hessian_values[1] = 2. * objective_multiplier - 2. * multipliers[0];
         hessian_values[2] = objective_multiplier;
         hessian_values[3] = 2. * objective_multiplier;
         hessian_values[4] = 12. * objective_multiplier * std::pow(x[3] - 1., 2);
         hessian_values[5] = 2. * objective_multiplier;
      }
      void compute_jacobian_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_jacobian_transposed_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
            const Vector<double>& multipliers, View<double> result) const override {
         result[0] = 2. * objective_multiplier * vector[0] + objective_multiplier * vector[1];
         result[1] = objective_multiplier * vector[0] + (2. * objective_multiplier - 2. * multipliers[0]) * vector[1];
         result[2] = 2. * objective_multiplier * vector[2];
         result[3] = 12. * objective_multiplier * std::pow(x[3] - 1., 2) * vector[3];
         result[4] = 2. * objective_multiplier * vector[4];
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
         x[2] = -0.2;
         x[3] = 2.;
         x[4] = 1.;
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      size_t number_jacobian_nonzeros() const override { return 7; }
      size_t number_hessian_nonzeros() const override { return 6; }

      size_t number_model_objective_evaluations() const override { return 0; }
      size_t number_model_constraints_evaluations() const override { return 0; }
      size_t number_model_objective_gradient_evaluations() const override { return 0; }
      size_t number_model_jacobian_evaluations() const override { return 0; }
      size_t number_model_hessian_evaluations() const override { return 0; }
      void reset_number_evaluations() const override { }

   private:
      const std::vector<uno_int> jacobian_rows{0, 0, 0, 1, 1, 2, 2};
      const std::vector<uno_int> jacobian_columns{0, 1, 2, 0, 1, 2, 3};
      const std::vector<double> variables_lower_bounds{0., -5., -INFINITY, -INFINITY, 0.5};
      const std::vector<double> variables_upper_bounds{10., INFINITY, 4., INFINITY, 0.5};
      const std::vector<double> constraints_lower_bounds{1., -INFINITY, -1.};
      const std::vector<double> constraints_upper_bounds{1., 3., 2.};
      const std::vector<size_t> equality_constraints{0};
      const std::vector<size_t> inequality_constraints{1, 2};
      const IntegerRange linear_constraints{1, 3};
      const IntegerRange nonlinear_constraints{0, 1};
   };
   // the ingredients of a subproblem at a nontrivial iterate
   template <typename Problem, typename InertiaCorrection = PrimalDualInertiaCorrection>
   struct Setting {
      static constexpr double residual_scaling_threshold = 100.;
      BarrierTestModel model{};
      Options options;
      OptimizationProblem problem{model, residual_scaling_threshold};
      InteriorPointParameters parameters{0.99, 1e10, 0.25, 10., 1e-2, 1e-2, 1e-5, 1., 1e-12};
      double barrier_parameter{0.1};
      std::unique_ptr<OptimizationProblem> barrier_problem;
      ExactHessian hessian_model{model};
      InertiaCorrection inertia_correction;
      std::unique_ptr<Subproblem> subproblem;
      Statistics statistics{};
      Evaluations evaluations{model};
      Iterate iterate{};

      Setting(): options(make_options()), inertia_correction(this->options) {
         if constexpr (std::is_same_v<Problem, PrimalDualInteriorPointProblem>) {
            this->barrier_problem = std::make_unique<PrimalDualInteriorPointProblem>(this->problem, this->parameters,
               this->barrier_parameter, residual_scaling_threshold);
         }
         const OptimizationProblem& subproblem_problem = this->barrier_problem ? *this->barrier_problem : this->problem;
         this->subproblem = std::make_unique<Subproblem>(subproblem_problem, this->hessian_model, this->inertia_correction);
         this->hessian_model.initialize_statistics(this->statistics);
         this->inertia_correction.initialize_statistics(this->statistics);
         this->iterate = Iterate(this->model.number_variables, this->model.number_constraints);
         this->model.initial_primal_point(this->iterate.primals);
         subproblem_problem.create_initial_iterate(this->iterate, this->evaluations);
         // nontrivial multipliers
         for (size_t index: Range(this->iterate.multipliers.constraints.size())) {
            this->iterate.multipliers.constraints[index] = 0.5 - 0.3 * static_cast<double>(index);
         }
         for (size_t index: Range(this->iterate.multipliers.lower_bounds.size())) {
            if (this->iterate.multipliers.lower_bounds[index] != 0.) this->iterate.multipliers.lower_bounds[index] = 0.7 + 0.1 * static_cast<double>(index);
            if (this->iterate.multipliers.upper_bounds[index] != 0.) this->iterate.multipliers.upper_bounds[index] = -0.4 - 0.2 * static_cast<double>(index);
         }
      }

      static Options make_options() {
         Options options;
         DefaultOptions::load(options);
         Presets::set(options, "ipopt");
         return options;
      }
   };

   // (problem, inertia correction strategy): with and without dual regularization
   template <typename P, typename I>
   struct Configuration {
      using Problem = P;
      using InertiaCorrection = I;
   };
   template <typename C>
   class NewtonSystemTest: public testing::Test { };
   using Configurations = testing::Types<Configuration<OptimizationProblem, PrimalDualInertiaCorrection>,
      Configuration<PrimalDualInteriorPointProblem, PrimalDualInertiaCorrection>,
      Configuration<OptimizationProblem, PrimalInertiaCorrection>,
      Configuration<PrimalDualInteriorPointProblem, PrimalInertiaCorrection>>;
   TYPED_TEST_SUITE(NewtonSystemTest, Configurations);
} // namespace

// the symbolic Newton system is the augmented matrix of the subproblem: pattern, block layout and values (bitwise)
TYPED_TEST(NewtonSystemTest, AssembledMatrixIsAugmentedMatrix) {
   Setting<typename TypeParam::Problem, typename TypeParam::InertiaCorrection> setting;
   const std::unique_ptr<NewtonSystem> newton_system = setting.subproblem->create_newton_system();

   COOLinearSystem reference(Indexing::Fortran_indexing);
   reference.initialize_augmented_system(*setting.subproblem);
   COOLinearSystem symbolic(Indexing::Fortran_indexing);
   symbolic.initialize(*newton_system);
   ASSERT_EQ(symbolic.dimension, reference.dimension);
   ASSERT_EQ(symbolic.number_nonzeros, reference.number_nonzeros);
   EXPECT_EQ(symbolic.matrix_row_indices, reference.matrix_row_indices);
   EXPECT_EQ(symbolic.matrix_column_indices, reference.matrix_column_indices);

   // block layout
   auto reference_blocks = setting.subproblem->compute_block_augmented_matrix(reference);
   const auto symbolic_blocks = newton_system->block_augmented_matrix(reference.matrix_values.view());
   EXPECT_EQ(symbolic_blocks.hessian.data(), reference_blocks.hessian.data());
   EXPECT_EQ(symbolic_blocks.hessian.size(), reference_blocks.hessian.size());
   EXPECT_EQ(symbolic_blocks.primal_inertia_correction.data(), reference_blocks.primal_inertia_correction.data());
   EXPECT_EQ(symbolic_blocks.primal_inertia_correction.size(), reference_blocks.primal_inertia_correction.size());
   EXPECT_EQ(symbolic_blocks.jacobian.data(), reference_blocks.jacobian.data());
   EXPECT_EQ(symbolic_blocks.jacobian.size(), reference_blocks.jacobian.size());
   EXPECT_EQ(symbolic_blocks.dual_inertia_correction.data(), reference_blocks.dual_inertia_correction.data());
   EXPECT_EQ(symbolic_blocks.dual_inertia_correction.size(), reference_blocks.dual_inertia_correction.size());

   // values (the regularization blocks are not evaluated)
   setting.subproblem->evaluate_lagrangian_hessian(setting.statistics, setting.iterate, reference_blocks.hessian);
   setting.subproblem->evaluate_jacobian(setting.iterate, reference_blocks.jacobian, setting.evaluations);
   reference_blocks.primal_inertia_correction.fill(0.);
   reference_blocks.dual_inertia_correction.fill(0.);
   symbolic.matrix_values.fill(0.);
   newton_system->evaluate({setting.statistics, setting.iterate, setting.evaluations}, symbolic.matrix_values.view());
   for (size_t nonzero_index: Range(reference.number_nonzeros)) {
      EXPECT_EQ(symbolic.matrix_values[nonzero_index], reference.matrix_values[nonzero_index]) << "entry " << nonzero_index;
   }
}

namespace {
   // y += A x with A given by its lower-triangular COO entries (symmetric) or by its full COO entries
   void add_coo_product(const std::vector<uno_int>& rows, const std::vector<uno_int>& columns, const Vector<double>& values,
         bool symmetric, const std::vector<double>& x, std::vector<double>& y) {
      for (size_t k: Range(values.size())) {
         const size_t i = static_cast<size_t>(rows[k]), j = static_cast<size_t>(columns[k]);
         y[i] += values[k] * x[j];
         if (symmetric && i != j) y[j] += values[k] * x[i];
      }
   }

   template <typename Symbol>
   void expect_operator_matches_entries(const Symbol& symbol, const NewtonSystemContext& context, size_t number_rows,
         size_t number_columns, bool symmetric) {
      const size_t number_nonzeros = uno::number_nonzeros(symbol);
      std::vector<uno_int> rows(number_nonzeros), columns(number_nonzeros);
      compute_sparsity(symbol, view(rows), view(columns), 0, 0, Indexing::C_indexing);
      Vector<double> values(number_nonzeros);
      evaluate(symbol, context, values.view());
      std::vector<double> x(number_columns);
      for (size_t index: Range(number_columns)) x[index] = 1. + 0.37 * static_cast<double>(index) - 0.05 * static_cast<double>(index * index);
      std::vector<double> reference(number_rows, 0.), result(number_rows, 0.);
      add_coo_product(rows, columns, values, symmetric, x, reference);
      add_matrix_vector_product(symbol, context, View<const double>(x.data(), x.size()), View<double>(result.data(), result.size()));
      for (size_t index: Range(number_rows)) {
         EXPECT_NEAR(result[index], reference[index], 1e-14 * (1. + std::abs(reference[index]))) << "row " << index;
      }
   }
} // namespace

// the operator products of the symbols coincide with the products with their assembled entries
TEST(NewtonSystemSymbols, OperatorProductsMatchEntries) {
   Setting<OptimizationProblem> setting;
   const NewtonSystemContext context{setting.statistics, setting.iterate, setting.evaluations};
   const size_t n = setting.problem.number_variables, m = setting.problem.number_constraints;
   setting.evaluations.evaluate_jacobian(setting.model, setting.iterate.primals);

   const LagrangianHessian hessian(setting.problem, setting.hessian_model, n);
   expect_operator_matches_entries(hessian, context, n, n, true);
   const ConstraintJacobian jacobian(setting.problem);
   expect_operator_matches_entries(jacobian, context, m, n, false);
   const std::vector<size_t> indices{0, 2, 3};
   const auto diagonal = diagonal_matrix(n, indices, [](size_t index) { return 0.5 + static_cast<double>(index); });
   expect_operator_matches_entries(diagonal, context, n, n, true);
   // sum: entries of both operands
   expect_operator_matches_entries(LagrangianHessian(setting.problem, setting.hessian_model, n) +
      diagonal_matrix(n, indices, [](size_t index) { return -1.5 * static_cast<double>(index); }), context, n, n, true);

   // transposed Jacobian product
   std::vector<double> v(m), reference(n, 0.), result(n, 0.);
   for (size_t index: Range(m)) v[index] = 0.3 - static_cast<double>(index);
   const size_t number_nonzeros = jacobian.number_nonzeros();
   std::vector<uno_int> rows(number_nonzeros), columns(number_nonzeros);
   jacobian.compute_sparsity(view(rows), view(columns), 0, 0, Indexing::C_indexing);
   Vector<double> values(number_nonzeros);
   jacobian.evaluate(context, values.view());
   for (size_t k: Range(number_nonzeros)) {
      reference[static_cast<size_t>(columns[k])] += values[k] * v[static_cast<size_t>(rows[k])];
   }
   jacobian.add_transposed_matrix_vector_product(context, View<const double>(v.data(), m), View<double>(result.data(), n));
   for (size_t index: Range(n)) {
      EXPECT_DOUBLE_EQ(result[index], reference[index]);
   }
}

// the structure of a symmetric block matrix: [A11 *; A21 A22] with offsets
TEST(NewtonSystemSymbols, SymmetricBlockMatrixLayout) {
   Setting<OptimizationProblem> setting;
   const size_t n = setting.problem.number_variables, m = setting.problem.number_constraints;
   const auto matrix = symmetric_block_matrix(Regularization(n), ConstraintJacobian(setting.problem), Regularization(m));
   EXPECT_EQ(number_rows(matrix), n + m);
   EXPECT_EQ(number_nonzeros(matrix), n + setting.problem.number_jacobian_nonzeros() + m);
   std::vector<uno_int> rows(number_nonzeros(matrix)), columns(number_nonzeros(matrix));
   compute_sparsity(matrix, view(rows), view(columns), 0, 0, Indexing::Fortran_indexing);
   // first entry of each block
   EXPECT_EQ(rows[0], 1);
   EXPECT_EQ(columns[0], 1);
   const size_t jacobian_start = n;
   EXPECT_EQ(rows[jacobian_start], static_cast<uno_int>(n) + setting.problem.get_jacobian_row_indices()[0] + 1);
   EXPECT_EQ(columns[jacobian_start], setting.problem.get_jacobian_column_indices()[0] + 1);
   const size_t dual_start = n + setting.problem.number_jacobian_nonzeros();
   EXPECT_EQ(rows[dual_start], static_cast<uno_int>(n) + 1);
   EXPECT_EQ(columns[dual_start], static_cast<uno_int>(n) + 1);
}
