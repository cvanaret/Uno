// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Symbolic Newton system: the matrices assembled from the symbols coincide (pattern and values, bitwise) with the
// augmented matrices of the subproblem (unconstrained/equality-constrained problem and primal-dual interior-point problem),
// and the operator products of the symbols coincide with the products with their assembled entries

#include <cmath>
#include <memory>
#include <stdexcept>
#include <type_traits>
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

namespace {
   // solve a dense linear system A x = b (Gaussian elimination with partial pivoting)
   std::vector<double> dense_solve(std::vector<std::vector<double>> A, std::vector<double> b) {
      const size_t n = b.size();
      for (size_t k: Range(n)) {
         size_t pivot = k;
         for (size_t i = k + 1; i < n; ++i) if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
         std::swap(A[k], A[pivot]);
         std::swap(b[k], b[pivot]);
         for (size_t i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / A[k][k];
            for (size_t j = k; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
         }
      }
      std::vector<double> x(n);
      for (size_t k = n; k-- > 0;) {
         double sum = b[k];
         for (size_t j = k + 1; j < n; ++j) sum -= A[k][j] * x[j];
         x[k] = sum / A[k][k];
      }
      return x;
   }
} // namespace

// condensation of the bound duals: the solution of the condensed system, completed by back substitution, solves the
// uncondensed Newton system of the barrier problem (stationarity, linearized constraints and linearized complementarity)
TEST(NewtonSystemCondensation, CondensedSolutionSolvesUncondensedSystem) {
   Setting<PrimalDualInteriorPointProblem> setting;
   const OptimizationProblem& barrier_problem = *setting.barrier_problem;
   const size_t n = barrier_problem.number_variables, m = barrier_problem.number_constraints;
   const size_t number_original_variables = setting.model.number_variables;
   ASSERT_EQ(n, size_t{7}); // 5 variables + 2 slacks
   ASSERT_EQ(m, size_t{4}); // 3 constraints + 1 fixed variable

   // internal bounds of the barrier problem: original variables (x4 is fixed: moved to the constraints), then the slacks
   const std::vector<double> lower_bounds{0., -5., -INFINITY, -INFINITY, -INFINITY, -INFINITY, -1.};
   const std::vector<double> upper_bounds{10., INFINITY, 4., INFINITY, INFINITY, 3., 2.};
   std::vector<size_t> lower_bounded, upper_bounded;
   Iterate& iterate = setting.iterate;
   for (size_t j: Range(n)) {
      if (std::isfinite(lower_bounds[j])) lower_bounded.push_back(j);
      if (std::isfinite(upper_bounds[j])) upper_bounded.push_back(j);
      iterate.multipliers.lower_bounds[j] = std::isfinite(lower_bounds[j]) ? 0.7 + 0.1 * static_cast<double>(j) : 0.;
      iterate.multipliers.upper_bounds[j] = std::isfinite(upper_bounds[j]) ? -0.4 - 0.2 * static_cast<double>(j) : 0.;
   }
   const size_t number_lower = lower_bounded.size(), number_upper = upper_bounded.size();

   // condensed system
   const std::unique_ptr<NewtonSystem> newton_system = setting.subproblem->create_newton_system();
   ASSERT_EQ(newton_system->dimension(), n + m);
   ASSERT_EQ(newton_system->number_eliminated(), number_lower + number_upper);
   const NewtonSystemContext context{setting.statistics, iterate, setting.evaluations};
   const size_t number_nonzeros = newton_system->number_nonzeros();
   std::vector<uno_int> rows(number_nonzeros), columns(number_nonzeros);
   newton_system->compute_sparsity(view(rows), view(columns), Indexing::C_indexing);
   Vector<double> values(number_nonzeros, 0.); // regularization blocks: 0
   newton_system->evaluate(context, values.view());
   std::vector<std::vector<double>> matrix(n + m, std::vector<double>(n + m, 0.));
   for (size_t k: Range(number_nonzeros)) {
      const size_t i = static_cast<size_t>(rows[k]), j = static_cast<size_t>(columns[k]);
      matrix[i][j] += values[k];
      if (i != j) matrix[j][i] += values[k];
   }
   Vector<double> rhs(n + m);
   newton_system->assemble_rhs(context, rhs.view());
   const std::vector<double> solution_vector = dense_solve(matrix, std::vector<double>(rhs.data(), rhs.data() + n + m));
   Vector<double> solution(n + m);
   for (size_t i: Range(n + m)) solution[i] = solution_vector[i];
   Vector<double> eliminated(number_lower + number_upper, 0.);
   newton_system->back_substitute(solution.view(), eliminated.view());

   // direction: Δx, Δy (the solution contains -Δy), Δz_L, Δz_U (scattered to the variables)
   const auto dx = view(solution.data(), n);
   std::vector<double> minus_dy(solution.data() + n, solution.data() + n + m);
   std::vector<double> dz_lower(n, 0.), dz_upper(n, 0.);
   for (size_t k: Range(number_lower)) dz_lower[lower_bounded[k]] = eliminated[k];
   for (size_t k: Range(number_upper)) dz_upper[upper_bounded[k]] = eliminated[number_lower + k];
   const double mu = setting.barrier_parameter;

   // independent check of the equations
   // ∇f + damping, recovered from the barrier gradient ∇f - μ/S_L - μ/S_U + damping
   Vector<double> gradient(n);
   barrier_problem.evaluate_objective_gradient(iterate, gradient.view(), setting.evaluations);
   for (size_t j: Range(n)) {
      if (std::isfinite(lower_bounds[j])) gradient[j] += mu / (iterate.primals[j] - lower_bounds[j]);
      if (std::isfinite(upper_bounds[j])) gradient[j] += mu / (iterate.primals[j] - upper_bounds[j]);
   }

   // stationarity: ∇²L Δx + J^T (-Δy) - Δz_L - Δz_U = -(∇f + damping - J^T y - z_L - z_U)
   Vector<double> hessian_product(number_original_variables, 0.);
   setting.problem.compute_hessian_vector_product(setting.hessian_model, view(iterate.primals.data(), number_original_variables),
      view(dx.data(), number_original_variables), iterate.multipliers, hessian_product.view());
   std::vector<double> stationarity(n, 0.), jacobian_y(n, 0.);
   barrier_problem.add_jacobian_transposed_vector_product(View<const double>(minus_dy.data(), m), View<double>(stationarity.data(), n),
      setting.evaluations);
   barrier_problem.add_jacobian_transposed_vector_product(iterate.multipliers.constraints.view(), View<double>(jacobian_y.data(), n),
      setting.evaluations);
   for (size_t j: Range(n)) {
      if (j < number_original_variables) stationarity[j] += hessian_product[j];
      stationarity[j] += -dz_lower[j] - dz_upper[j];
      const double expected = -(gradient[j] - jacobian_y[j] - iterate.multipliers.lower_bounds[j] - iterate.multipliers.upper_bounds[j]);
      EXPECT_NEAR(stationarity[j], expected, 1e-12 * (1. + std::abs(expected))) << "stationarity " << j;
   }

   // linearized constraints: J Δx = -(c(x) - c_L)
   std::vector<double> constraints(m, 0.), jacobian_product(m, 0.);
   barrier_problem.evaluate_constraints(iterate, View<double>(constraints.data(), m), setting.evaluations);
   barrier_problem.add_jacobian_vector_product(View<const double>(dx.data(), n), View<double>(jacobian_product.data(), m), setting.evaluations);
   for (size_t i: Range(m)) {
      const double expected = -(constraints[i] - barrier_problem.get_constraints_lower_bounds()[i]);
      EXPECT_NEAR(jacobian_product[i], expected, 1e-12 * (1. + std::abs(expected))) << "constraint " << i;
   }

   // linearized complementarity: z Δx + S Δz = μ - S z
   for (const size_t j: lower_bounded) {
      const double distance = iterate.primals[j] - lower_bounds[j], z = iterate.multipliers.lower_bounds[j];
      EXPECT_NEAR(z * dx[j] + distance * dz_lower[j], mu - distance * z, 1e-12) << "lower complementarity " << j;
   }
   for (const size_t j: upper_bounded) {
      const double distance = iterate.primals[j] - upper_bounds[j], z = iterate.multipliers.upper_bounds[j];
      EXPECT_NEAR(z * dx[j] + distance * dz_upper[j], mu - distance * z, 1e-12) << "upper complementarity " << j;
   }

   // the problem recovers the same bound duals from the primal direction (linearized complementarity)
   Multipliers direction_multipliers(n, m);
   double step_length = 0.;
   Vector<double> primal_direction(n);
   for (size_t j: Range(n)) primal_direction[j] = dx[j];
   barrier_problem.compute_bound_dual_direction(iterate.primals, iterate.multipliers, primal_direction, direction_multipliers, step_length);
   for (size_t j: Range(n)) {
      EXPECT_EQ(direction_multipliers.lower_bounds[j], dz_lower[j]);
      EXPECT_EQ(direction_multipliers.upper_bounds[j], dz_upper[j]);
   }
}

// generic condensation of a small block system: kept unknowns u (block 0, dimension 2) and w (block 1, dimension 1),
// eliminated unknowns v (block 2, dimension 2) coupled with u0 and u1, and t (block 3, dimension 1) coupled with u1
TEST(NewtonSystemCondensation, SmallBlockSystem) {
   Setting<OptimizationProblem> setting;
   const NewtonSystemContext context{setting.statistics, setting.iterate, setting.evaluations};
   const auto system = kkt_system(
      block_matrix({2, 1, 2, 1},
         block<0, 0>(diagonal(2, [](size_t i) { return 3. + static_cast<double>(i); })),
         block<0, 1>(scaled_selection(2, 1, 1, ConstantIndex{1}, IdentityIndex{}, ConstantValue{2.})),
         block<1, 0>(scaled_selection(1, 2, 1, IdentityIndex{}, ConstantIndex{1}, ConstantValue{2.})),
         block<0, 2>(scaled_selection(2, 2, 2, IdentityIndex{}, IdentityIndex{}, [](size_t i) { return i == 0 ? -1. : 0.5; })),
         block<2, 0>(scaled_selection(2, 2, 2, IdentityIndex{}, IdentityIndex{}, [](size_t i) { return i == 0 ? 2. : -1.; })),
         block<2, 2>(diagonal(2, [](size_t i) { return i == 0 ? 4. : -2.; })),
         block<0, 3>(scaled_selection(2, 1, 1, ConstantIndex{1}, IdentityIndex{}, ConstantValue{-3.})),
         block<3, 0>(scaled_selection(1, 2, 1, IdentityIndex{}, ConstantIndex{1}, ConstantValue{1.5})),
         block<3, 3>(diagonal(1, ConstantValue{5.}))),
      block_vector(component_vector(2, [](size_t i) { return 1. - static_cast<double>(i); }), component_vector(1, ConstantValue{0.25}),
         component_vector(2, [](size_t i) { return i == 0 ? 0.5 : 3.; }), component_vector(1, ConstantValue{-1.})));
   Condensation condensation(system, {2, 3});
   ASSERT_EQ(condensation.number_kept(), size_t{3});
   ASSERT_EQ(condensation.number_eliminated(), size_t{3});
   ASSERT_EQ(condensation.number_diagonal_entries(), size_t{2});
   condensation.evaluate_matrix(system);
   std::vector<double> rhs(3);
   condensation.evaluate_rhs(system, context, View<double>(rhs.data(), 3));

   // condensed system (dense): kept blocks + diagonal Schur complement
   std::vector<std::vector<double>> matrix{{3., 0., 0.}, {0., 4., 2.}, {0., 2., 0.}};
   for (size_t position: Range(condensation.number_diagonal_entries())) {
      const size_t index = condensation.diagonal_index(position);
      matrix[index][index] += condensation.diagonal_value(position);
   }
   const std::vector<double> kept = dense_solve(matrix, rhs);
   std::vector<double> eliminated(3);
   condensation.back_substitute(View<const double>(kept.data(), 3), View<double>(eliminated.data(), 3));

   // the full solution solves the uncondensed system
   std::vector<double> full{kept[0], kept[1], kept[2], eliminated[0], eliminated[1], eliminated[2]}, product(6, 0.), full_rhs(6);
   system.get_matrix().add_matrix_vector_product(context, View<const double>(full.data(), 6), View<double>(product.data(), 6));
   system.evaluate_rhs(context, View<double>(full_rhs.data(), 6));
   for (size_t i: Range(6)) {
      EXPECT_NEAR(product[i], full_rhs[i], 1e-14) << "row " << i;
   }
}

// blocks that cannot be condensed into a diagonal Schur complement are rejected
TEST(NewtonSystemCondensation, RejectsNonDiagonalSchurComplement) {
   // v0 coupled with u0 in the kept row, but with u1 in its own row
   const auto system = kkt_system(
      block_matrix({2, 0, 1},
         block<0, 2>(scaled_selection(2, 1, 1, ConstantIndex{0}, IdentityIndex{}, ConstantValue{1.})),
         block<2, 0>(scaled_selection(1, 2, 1, IdentityIndex{}, ConstantIndex{1}, ConstantValue{1.})),
         block<2, 2>(diagonal(1, ConstantValue{1.}))),
      block_vector(component_vector(2, ConstantValue{0.}), component_vector(0, ConstantValue{0.}), component_vector(1, ConstantValue{0.})));
   EXPECT_THROW(Condensation(system, {2}), std::invalid_argument);
}

// symbolic cancellation: the right-hand side in residual form (linear terms in the current values of the eliminated
// unknowns) yields the same direction as the plain right-hand side, without evaluating the cancelled terms
TEST(NewtonSystemCondensation, ResidualFormCancelsSymbolically) {
   Setting<OptimizationProblem> setting;
   const NewtonSystemContext context{setting.statistics, setting.iterate, setting.evaluations};
   // kept u (block 0, dimension 2), dual w (block 1, dimension 1), eliminated v (block 2, dimension 2) with current value v0
   const auto v0 = [](size_t i) { return i == 0 ? 1e8 : -3.; }; // large current value: catastrophic cancellation if evaluated
   const auto matrix = [] {
      return block_matrix({2, 1, 2},
         block<0, 0>(diagonal(2, [](size_t i) { return 3. + static_cast<double>(i); })),
         block<0, 1>(scaled_selection(2, 1, 1, ConstantIndex{1}, IdentityIndex{}, ConstantValue{2.})),
         block<1, 0>(scaled_selection(1, 2, 1, IdentityIndex{}, ConstantIndex{1}, ConstantValue{2.})),
         block<0, 2>(scaled_selection(2, 2, 2, IdentityIndex{}, IdentityIndex{}, ConstantValue{-1.})),
         block<2, 0>(scaled_selection(2, 2, 2, IdentityIndex{}, IdentityIndex{}, [](size_t i) { return i == 0 ? 2. : 0.5; })),
         block<2, 2>(diagonal(2, [](size_t i) { return i == 0 ? 1e-3 : 4.; })));
   };
   // residual form: kept r = c_u - K(0,2) v0, eliminated r = q - K(2,2) v0
   const auto c_u = [](size_t i) { return i == 0 ? 0.1 : -1.; };
   const auto q = [](size_t i) { return i == 0 ? 0.3 : 2.; };
   const auto residual_system = kkt_system(matrix(), block_vector(newton_rhs(component_vector(2, c_u), linear_term<2>(v0)),
      component_vector(1, ConstantValue{0.75}), newton_rhs(component_vector(2, q), linear_term<2>(v0))));
   // plain form: the same right-hand side, evaluated
   const auto plain_system = kkt_system(matrix(), block_vector(
      component_vector(2, [&](size_t i) { return c_u(i) + v0(i); }), component_vector(1, ConstantValue{0.75}),
      component_vector(2, [&](size_t i) { return q(i) - (i == 0 ? 1e-3 : 4.) * v0(i); })));

   const auto solve = [&](const auto& system, std::vector<double>& eliminated) {
      Condensation condensation(system, {2});
      condensation.evaluate_matrix(system);
      std::vector<double> rhs(3);
      condensation.evaluate_rhs(system, context, View<double>(rhs.data(), 3));
      std::vector<std::vector<double>> dense{{3., 0., 0.}, {0., 4., 2.}, {0., 2., 0.}};
      for (size_t position: Range(condensation.number_diagonal_entries())) {
         const size_t index = condensation.diagonal_index(position);
         dense[index][index] += condensation.diagonal_value(position);
      }
      const std::vector<double> kept = dense_solve(dense, rhs);
      eliminated.resize(2);
      condensation.back_substitute(View<const double>(kept.data(), 3), View<double>(eliminated.data(), 2));
      return kept;
   };
   std::vector<double> residual_eliminated, plain_eliminated;
   const std::vector<double> residual_kept = solve(residual_system, residual_eliminated);
   const std::vector<double> plain_kept = solve(plain_system, plain_eliminated);

   // exact solution: the cancelled terms do not appear. Kept rows: u0 (no coupling with w), u1 and w
   // condensed: (3 + 2/1e-3) u0 = 0.1 + 0.3/1e-3
   const double u0 = (0.1 + 0.3 / 1e-3) / (3. + 2. / 1e-3);
   EXPECT_NEAR(residual_kept[0], u0, 1e-15 * std::abs(u0));
   EXPECT_NEAR(residual_eliminated[0], (0.3 - 2. * u0) / 1e-3 - 1e8, 1e-7);
   // the plain form loses accuracy in u0: the current value 1e8 is added in the kept row and subtracted by the correction
   EXPECT_GT(std::abs(plain_kept[0] - u0), 1e3 * std::abs(residual_kept[0] - u0) + 1e-14);
   // both forms agree on the well-scaled entries
   EXPECT_NEAR(residual_kept[1], plain_kept[1], 1e-12);
   EXPECT_NEAR(residual_kept[2], plain_kept[2], 1e-12);

   // the uncondensed right-hand side of the residual form contains the linear terms
   std::vector<double> full_rhs(5), plain_full_rhs(5);
   residual_system.evaluate_rhs(context, View<double>(full_rhs.data(), 5));
   plain_system.evaluate_rhs(context, View<double>(plain_full_rhs.data(), 5));
   for (size_t i: Range(5)) {
      EXPECT_DOUBLE_EQ(full_rhs[i], plain_full_rhs[i]) << "row " << i;
   }
}

// a kept row coupled with a condensed block in residual form must declare the linear term of the coupling
TEST(NewtonSystemCondensation, RejectsMissingLinearTerm) {
   const auto system = kkt_system(
      block_matrix({1, 0, 1},
         block<0, 2>(scaled_selection(1, 1, 1, IdentityIndex{}, IdentityIndex{}, ConstantValue{-1.})),
         block<2, 0>(scaled_selection(1, 1, 1, IdentityIndex{}, IdentityIndex{}, ConstantValue{1.})),
         block<2, 2>(diagonal(1, ConstantValue{1.}))),
      block_vector(component_vector(1, ConstantValue{0.}), component_vector(0, ConstantValue{0.}),
         newton_rhs(component_vector(1, ConstantValue{0.}), linear_term<2>(ConstantValue{1.}))));
   EXPECT_THROW(Condensation(system, {2}), std::invalid_argument);
}

// the dense declaration (rows with explicit zero blocks) is lowered to the sparse declaration: same block types, in row order
TEST(BlockMatrixDeclaration, RowsLowerToSparseBlocks) {
   const auto D = [] { return diagonal(2, ConstantValue{2.}); };
   const auto C = [] { return scaled_selection(1, 2, 1, IdentityIndex{}, ConstantIndex{1}, ConstantValue{3.}); };
   const auto Ct = [] { return scaled_selection(2, 1, 1, ConstantIndex{1}, IdentityIndex{}, ConstantValue{3.}); };
   const auto E = [] { return diagonal(1, ConstantValue{-1.}); };
   // clang-format off
   auto dense = block_matrix({2, 1, 1},
      row(D(),  Ct(),  zero),
      row(C(),  zero,  zero),
      row(zero, zero,  E()));
   // clang-format on
   auto sparse = block_matrix({2, 1, 1}, block<0, 0>(D()), block<0, 1>(Ct()), block<1, 0>(C()), block<2, 2>(E()));
   static_assert(std::is_same_v<decltype(dense), decltype(sparse)>);

   // same products
   Setting<OptimizationProblem> setting;
   const NewtonSystemContext context{setting.statistics, setting.iterate, setting.evaluations};
   const std::vector<double> x{1., -2., 0.5, 4.};
   std::vector<double> y_dense(4, 0.), y_sparse(4, 0.);
   dense.add_matrix_vector_product(context, View<const double>(x.data(), 4), View<double>(y_dense.data(), 4));
   sparse.add_matrix_vector_product(context, View<const double>(x.data(), 4), View<double>(y_sparse.data(), 4));
   EXPECT_EQ(y_dense, y_sparse);
   const std::vector<double> expected{2., -4. + 1.5, -6., -4.};
   EXPECT_EQ(y_dense, expected);

   // an lvalue entry is stored by reference, an rvalue entry is owned
   const auto lvalue_diagonal = D();
   auto with_reference = block_matrix({2, 1}, row(lvalue_diagonal, zero), row(zero, E()));
   static_assert(std::is_same_v<decltype(with_reference),
      BlockMatrix<Block<0, 0, const std::decay_t<decltype(lvalue_diagonal)>&>, Block<1, 1, std::decay_t<decltype(E())>>>>);
   const auto* stored_diagonal = &with_reference.get_block<0, 0>();
   EXPECT_EQ(stored_diagonal, &lvalue_diagonal);
}

// the number of block dimensions must match the number of rows
TEST(BlockMatrixDeclaration, RejectsDimensionMismatch) {
   EXPECT_THROW(block_matrix({2, 1, 1}, row(diagonal(2, ConstantValue{1.}), zero), row(zero, diagonal(1, ConstantValue{1.}))),
      std::invalid_argument);
}

// back_substitute_block only accepts the own pivot term as linear term
TEST(NewtonSystemCondensation, OwnPivotTermTrait) {
   using Constant = ComponentVector<ConstantValue>;
   static_assert(has_only_own_pivot_term_v<2, Constant>);
   static_assert(has_only_own_pivot_term_v<2, NewtonRHS<Constant>>);
   static_assert(has_only_own_pivot_term_v<2, NewtonRHS<Constant, LinearTerm<2, ConstantValue>>>);
   static_assert(!has_only_own_pivot_term_v<2, NewtonRHS<Constant, LinearTerm<3, ConstantValue>>>);
   static_assert(!has_only_own_pivot_term_v<2, NewtonRHS<Constant, LinearTerm<2, ConstantValue>, LinearTerm<2, ConstantValue>>>);
}

#ifndef NDEBUG
// debug check: the linear term of a kept row cancelled symbolically must be the same vector as the own pivot term
TEST(NewtonSystemCondensation, DetectsInconsistentCancelledTerms) {
   Setting<OptimizationProblem> setting;
   const NewtonSystemContext context{setting.statistics, setting.iterate, setting.evaluations};
   const auto make_system = [](double kept_row_value) {
      // clang-format off
      return kkt_system(
         block_matrix({1, 0, 1},
            row(diagonal(1, ConstantValue{1.}),                                               zero, scaled_selection(1, 1, 1, IdentityIndex{}, IdentityIndex{}, ConstantValue{-1.})),
            row(zero,                                                                         zero, zero),
            row(scaled_selection(1, 1, 1, IdentityIndex{}, IdentityIndex{}, ConstantValue{2.}), zero, diagonal(1, ConstantValue{4.}))),
         block_vector(newton_rhs(component_vector(1, ConstantValue{0.}), linear_term<2>(ConstantValue{kept_row_value})),
            component_vector(0, ConstantValue{0.}),
            newton_rhs(component_vector(1, ConstantValue{1.}), linear_term<2>(ConstantValue{0.5}))));
      // clang-format on
   };
   std::vector<double> rhs(1);
   const auto consistent = make_system(0.5);
   Condensation condensation(consistent, {2});
   condensation.evaluate_matrix(consistent);
   EXPECT_NO_THROW(condensation.evaluate_rhs(consistent, context, View<double>(rhs.data(), 1)));

   const auto inconsistent = make_system(0.25);
   Condensation inconsistent_condensation(inconsistent, {2});
   inconsistent_condensation.evaluate_matrix(inconsistent);
   EXPECT_THROW(inconsistent_condensation.evaluate_rhs(inconsistent, context, View<double>(rhs.data(), 1)), std::logic_error);
}
#endif
