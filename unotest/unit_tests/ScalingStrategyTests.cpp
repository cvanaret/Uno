// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// linear system scaling strategy, applied by DirectSymmetricIndefiniteLinearSolver: on-demand activation, at most one scaling
// computation per matrix, consistency between factorization and solves. The solver is a dense mock that records the values
// it factorizes

#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/COOLinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/ScalingStrategy.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#if defined(HAS_HSL) || defined(HAS_MA57)
#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#endif
#ifdef HAS_MUMPS
#include "ingredients/subproblem_solvers/MUMPS/MUMPSSolver.hpp"
#endif
#include <algorithm>
#include <string>
#include "linear_algebra/Indexing.hpp"

using namespace uno;

namespace {
   // dense LU with partial pivoting of the symmetric COO matrix
   class MockSolver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      COOLinearSystem linear_system{Indexing::C_indexing};
      std::vector<double> factorized_values{};
      size_t number_factorizations{0};
      std::optional<size_t> delayed_pivots{0};
      bool throw_in_factorization{false};

      void initialize_memory() override { }
      void do_symbolic_analysis() override { }
      [[nodiscard]] Inertia get_inertia() const override { return {this->linear_system.dimension, 0, 0}; }
      [[nodiscard]] size_t number_negative_eigenvalues() const override { return 0; }
      [[nodiscard]] bool matrix_is_singular() const override { return false; }
      [[nodiscard]] size_t rank() const override { return this->linear_system.dimension; }
      [[nodiscard]] std::optional<size_t> number_delayed_pivots() const override { return this->delayed_pivots; }
      [[nodiscard]] LinearSystem& get_linear_system() override { return this->linear_system; }
      [[nodiscard]] COOLinearSystem& get_coo_linear_system() override { return this->linear_system; }
      [[nodiscard]] ScalingStrategy& strategy() { return *this->scaling_strategy; }

   protected:
      void factorize(bool /*is_matrix_positive_definite*/) override {
         ++this->number_factorizations;
         const size_t n = this->linear_system.dimension;
         const double* values = this->linear_system.matrix_values.data();
         this->factorized_values.assign(values, values + this->linear_system.number_nonzeros);
         if (this->throw_in_factorization) {
            throw std::runtime_error("factorization failure");
         }
         this->lu.assign(n * n, 0.);
         for (size_t k = 0; k < this->linear_system.number_nonzeros; ++k) {
            const auto i = static_cast<size_t>(this->linear_system.matrix_row_indices[k]);
            const auto j = static_cast<size_t>(this->linear_system.matrix_column_indices[k]);
            this->lu[i * n + j] += values[k];
            if (i != j) {
               this->lu[j * n + i] += values[k];
            }
         }
         this->pivots.resize(n);
         for (size_t k = 0; k < n; ++k) {
            size_t pivot = k;
            for (size_t i = k + 1; i < n; ++i) {
               if (std::abs(this->lu[i * n + k]) > std::abs(this->lu[pivot * n + k])) {
                  pivot = i;
               }
            }
            this->pivots[k] = pivot;
            for (size_t j = 0; j < n; ++j) {
               std::swap(this->lu[k * n + j], this->lu[pivot * n + j]);
            }
            for (size_t i = k + 1; i < n; ++i) {
               this->lu[i * n + k] /= this->lu[k * n + k];
               for (size_t j = k + 1; j < n; ++j) {
                  this->lu[i * n + j] -= this->lu[i * n + k] * this->lu[k * n + j];
               }
            }
         }
      }
      void solve(const double* rhs, double* solution, size_t number_of_rhs) override {
         const size_t n = this->linear_system.dimension;
         std::vector<double> block(rhs, rhs + n * number_of_rhs); // rhs and solution may alias
         for (size_t column = 0; column < number_of_rhs; ++column) {
            double* x = block.data() + column * n;
            for (size_t k = 0; k < n; ++k) {
               std::swap(x[k], x[this->pivots[k]]);
               for (size_t i = k + 1; i < n; ++i) {
                  x[i] -= this->lu[i * n + k] * x[k];
               }
            }
            for (size_t k = n; k-- > 0;) {
               for (size_t j = k + 1; j < n; ++j) {
                  x[k] -= this->lu[k * n + j] * x[j];
               }
               x[k] /= this->lu[k * n + k];
            }
         }
         std::copy(block.begin(), block.end(), solution);
      }

   private:
      std::vector<double> lu{};
      std::vector<size_t> pivots{};
   };

   // badly scaled 4x4 KKT matrix [H J^T; J 0] (lower triangle), the last entry is a primal regularization of x0
   void fill(COOLinearSystem& linear_system) {
      const std::vector<uno_int> rows{0, 1, 1, 2, 2, 3, 0};
      const std::vector<uno_int> columns{0, 0, 1, 0, 1, 1, 0};
      const std::vector<double> values{1e6, 1., 1e-4, 1e3, 2e-3, 5., 0.};
      linear_system.dimension = 4;
      linear_system.number_nonzeros = values.size();
      linear_system.matrix_row_indices = rows;
      linear_system.matrix_column_indices = columns;
      linear_system.matrix_values = Vector<double>(values.size());
      for (size_t index = 0; index < values.size(); ++index) {
         linear_system.matrix_values[index] = values[index];
      }
      linear_system.rhs = Vector<double>{1., -2., 3., 4.};
      linear_system.solution = Vector<double>(4);
   }

   struct Setting {
      MockSolver solver{};
      MockSolver* mock{&this->solver};

      explicit Setting(const LinearSystemScalingParameters& parameters) {
         this->solver.set_scaling_strategy(std::make_unique<ScalingStrategy>(parameters));
         fill(this->solver.linear_system);
         this->solver.initialize_memory();
         this->solver.do_symbolic_analysis();
      }

      [[nodiscard]] std::vector<double> values() const {
         const auto& values = this->solver.linear_system.matrix_values;
         return {values.data(), values.data() + values.size()};
      }
      [[nodiscard]] ScalingStrategy& strategy() { return this->solver.strategy(); }

      // first factorization of a new matrix
      void factorize_new_matrix() {
         this->solver.begin_new_matrix();
         this->solver.do_numerical_factorization(false);
      }
   };

   LinearSystemScalingParameters on_demand(LinearSystemScalingMethod method) {
      return {method, true, 0.05};
   }

   // solution of the unscaled system K x = b with the unscaled mock
   std::vector<double> reference_solution() {
      MockSolver mock;
      fill(mock.linear_system);
      mock.do_numerical_factorization(false);
      std::vector<double> solution(4);
      mock.solve_indefinite_system(solution.data());
      return solution;
   }

   void expect_near(const std::vector<double>& x, const std::vector<double>& y) {
      ASSERT_EQ(x.size(), y.size());
      for (size_t index = 0; index < x.size(); ++index) {
         EXPECT_NEAR(x[index], y[index], 1e-9 * (1. + std::abs(y[index])));
      }
   }
} // namespace

TEST(ScalingStrategy, AlwaysOnWhenNotOnDemand) {
   Setting setting({LinearSystemScalingMethod::RUIZ, false, 0.05});
   EXPECT_TRUE(setting.strategy().is_active());
   setting.factorize_new_matrix();
   EXPECT_TRUE(setting.strategy().is_last_factorization_scaled());
   EXPECT_EQ(setting.strategy().get_number_scaling_computations(), 1);
   EXPECT_FALSE(setting.solver.activate_scaling()); // already on
}

TEST(ScalingStrategy, NoneCannotBeActivated) {
   Setting setting({LinearSystemScalingMethod::NONE, true, 0.});
   EXPECT_FALSE(setting.solver.activate_scaling());
   setting.mock->delayed_pivots = 4;
   setting.factorize_new_matrix();
   setting.factorize_new_matrix();
   EXPECT_FALSE(setting.strategy().is_active());
}

// once switched on, the factorized matrix is S K S, the callers still see K, and the solutions are those of K x = b
TEST(ScalingStrategy, ScaledFactorizationAndSolve) {
   Setting setting(on_demand(LinearSystemScalingMethod::RUIZ));
   const std::vector<double> original = setting.values();
   setting.factorize_new_matrix();
   ASSERT_TRUE(setting.solver.activate_scaling());
   EXPECT_FALSE(setting.solver.activate_scaling()); // sticky
   EXPECT_TRUE(setting.strategy().is_active());
   // not refactorized yet: the solves are consistent with the (unscaled) factorization
   std::vector<double> solution(4);
   setting.solver.solve_indefinite_system(solution.data());
   expect_near(solution, reference_solution());

   setting.solver.do_numerical_factorization(false); // refactorization of the same matrix
   EXPECT_TRUE(setting.strategy().is_last_factorization_scaled());
   EXPECT_EQ(setting.strategy().get_number_scaling_computations(), 1);
   EXPECT_EQ(setting.values(), original); // restored bitwise
   const auto& scaling = setting.strategy().get_scaling();
   const auto& linear_system = setting.mock->linear_system;
   for (size_t index = 0; index < original.size(); ++index) {
      const double expected = scaling[static_cast<size_t>(linear_system.matrix_row_indices[index])] * original[index] *
         scaling[static_cast<size_t>(linear_system.matrix_column_indices[index])];
      EXPECT_DOUBLE_EQ(setting.mock->factorized_values[index], expected);
   }
   // single RHS (from the linear system)
   setting.solver.solve_indefinite_system(solution.data());
   expect_near(solution, reference_solution());
   // multiple RHS, aliased with the solution
   std::vector<double> block{1., -2., 3., 4., 2., -4., 6., 8.};
   setting.solver.solve_indefinite_system(block.data(), block.data(), 2);
   const std::vector<double> reference = reference_solution();
   for (size_t index = 0; index < 4; ++index) {
      EXPECT_NEAR(block[index], reference[index], 1e-9 * (1. + std::abs(reference[index])));
      EXPECT_NEAR(block[4 + index], 2. * reference[index], 1e-9 * (1. + std::abs(reference[index])));
   }
}

// without a strategy, the solver never scales and cannot be switched to scaling
TEST(ScalingStrategy, NoStrategy) {
   MockSolver solver;
   fill(solver.linear_system);
   const std::vector<double> original(solver.linear_system.matrix_values.data(),
      solver.linear_system.matrix_values.data() + solver.linear_system.number_nonzeros);
   EXPECT_EQ(solver.get_scaling_strategy(), nullptr);
   EXPECT_FALSE(solver.activate_scaling());
   solver.begin_new_matrix();
   solver.do_numerical_factorization(false);
   EXPECT_EQ(solver.factorized_values, original);
   std::vector<double> solution(4);
   solver.solve_indefinite_system(solution.data());
   expect_near(solution, reference_solution());
}

namespace {
   bool is_available(const std::string& linear_solver) {
      const auto solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      return std::find(solvers.begin(), solvers.end(), linear_solver) != solvers.end();
   }

   Options default_options(const std::string& linear_solver) {
      Options options;
      DefaultOptions::load(options);
      options.set_string("linear_solver", linear_solver);
      return options;
   }

   // whether the solver scales the matrix itself (MA57: MC64, MUMPS: ICNTL(8))
   std::optional<bool> uses_internal_scaling(const DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) {
#if defined(HAS_HSL) || defined(HAS_MA57)
      if (const auto* ma57 = dynamic_cast<const MA57Solver*>(&linear_solver)) {
         return ma57->uses_internal_scaling();
      }
#endif
#ifdef HAS_MUMPS
      if (const auto* mumps = dynamic_cast<const MUMPSSolver*>(&linear_solver)) {
         return mumps->uses_internal_scaling();
      }
#endif
      (void)linear_solver;
      return std::nullopt;
   }
} // namespace

// the solvers' own scalings (MA57: MC64, MUMPS: ICNTL(8)) follow MA57_use_scaling and MUMPS_use_scaling,
// independently of linear_system_scaling
TEST(ScalingStrategy, DefaultNativeScaling) {
   for (const std::string& linear_solver: {std::string("MA57"), std::string("MUMPS")}) {
      if (!is_available(linear_solver)) {
         continue;
      }
      const Options options = default_options(linear_solver);
      const bool expected = options.get_bool(linear_solver + "_use_scaling");
      EXPECT_EQ(uses_internal_scaling(*SymmetricIndefiniteLinearSolverFactory::create(options)), expected) << linear_solver;
   }
}
