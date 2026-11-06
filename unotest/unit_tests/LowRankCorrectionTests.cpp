// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <limits>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/hessian_models/quasi_newton/direct/DirectQuasiNewtonHessian.hpp"
#include "ingredients/inertia_correction_strategies/NoInertiaCorrection.hpp"
#include "ingredients/inertia_correction_strategies/PrimalDualInertiaCorrection.hpp"
#include "ingredients/subproblem_solvers/EQPSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Statistics.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/COOLinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/LowRankCorrection.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   class SyntheticHessian: public DirectQuasiNewtonHessian {
   public:
      SyntheticHessian(const Model& model, const Options& options): DirectQuasiNewtonHessian("synthetic", model, 1., options) { }
      std::vector<std::vector<double>> columns{};
      std::vector<double> coefficients{};
      size_t get_correction_rank() const override { return this->coefficients.size(); }
      View<const double> get_correction_column(size_t k) const override { return view(this->columns[k].data(), this->columns[k].size()); }
      double get_correction_column_scaling(size_t k) const override { return this->coefficients[k]; }
      bool is_positive_definite() const override { return false; }
      void initialize_statistics(Statistics&) const override { }
      void notify_trial_iterate(Statistics&, const Iterate&, const Iterate&, Evaluations&, Evaluations&) override { }
      void compute_hessian_vector_product(View<const double>, View<const double> v, double,
            const Vector<double>&, View<double> result) override {
         result = v;
         for (size_t k = 0; k < this->coefficients.size(); ++k) {
            double dot = 0.;
            for (size_t i = 0; i < v.size(); ++i) dot += this->columns[k][i] * v[i];
            for (size_t i = 0; i < v.size(); ++i) result[i] += this->coefficients[k] * this->columns[k][i] * dot;
         }
      }
   protected:
      void recompute_hessian_representation() override { }
   };

   class DiagonalSolver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      std::vector<double> diagonal{1., 1.};
      COOLinearSystem system{0};
      DiagonalSolver() { this->system.dimension = 2; this->system.rhs.resize(2); }
      void initialize_memory() override { }
      void do_symbolic_analysis() override { }
      void do_numerical_factorization(bool) override { }
      void solve_indefinite_system(double* solution) override { this->solve_indefinite_system(this->system.rhs.data(), solution, 1); }
      void solve_indefinite_system(const double* rhs, double* solution, size_t nrhs) override {
         for (size_t k = 0; k < nrhs * this->system.dimension; ++k) solution[k] = rhs[k] / this->diagonal[k % this->system.dimension];
      }
      Inertia get_inertia() const override {
         const size_t positive = static_cast<size_t>(std::count_if(this->diagonal.begin(), this->diagonal.end(), [](double d) { return d > 0.; }));
         return {positive, this->diagonal.size() - positive, 0};
      }
      size_t number_negative_eigenvalues() const override { return this->get_inertia().negative; }
      bool matrix_is_singular() const override { return false; }
      size_t rank() const override { return this->diagonal.size(); }
      LinearSystem& get_linear_system() override { return this->system; }
   };

   class LowRankCorrectionTest: public testing::Test {
   protected:
      Options options{};
      DenseTestModel model{"quadratic", 2, 0, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return (x[0] * x[0] + x[1] * x[1]) / 2.; },
         [](const Vector<double>& x, Vector<double>& g) { g = x; },
         [](const Vector<double>&, Vector<double>&) { }, [](const Vector<double>&, double*) { },
         [](const Vector<double>&, double w, const Vector<double>&, double* h) { h[0] = w; h[1] = 0.; h[2] = w; },
         {-Inf, -Inf}, {Inf, Inf}, {}, {}, {1., 2.}};
      OptimizationProblem problem{model, 100.};
      NoInertiaCorrection regularization{};
      DiagonalSolver solver{};
      std::unique_ptr<SyntheticHessian> hessian{};
      std::unique_ptr<Subproblem> subproblem{};
      std::unique_ptr<WoodburyCorrection> correction{};
      void SetUp() override {
         DefaultOptions::load(this->options);
         this->hessian = std::make_unique<SyntheticHessian>(this->model, this->options);
         this->subproblem = std::make_unique<Subproblem>(this->problem, *this->hessian, this->regularization);
         this->correction = std::make_unique<WoodburyCorrection>(*this->hessian);
      }
      Inertia update() {
         this->correction->update(*this->subproblem, this->solver);
         return this->correction->get_inertia(this->solver.get_inertia());
      }
   };
} // namespace

TEST_F(LowRankCorrectionTest, RankZeroPreservesBaseInertiaAndSolution) {
   EXPECT_TRUE(this->update() == Inertia(2, 0, 0));
   Vector<double> b{1., 2.};
   this->correction->apply(b);
   EXPECT_DOUBLE_EQ(b[0], 1.);
   EXPECT_DOUBLE_EQ(b[1], 2.);
}

TEST_F(LowRankCorrectionTest, NegativeRankOneTermChangesInertia) {
   this->hessian->columns = {{2., 0.}};
   this->hessian->coefficients = {-1.};
   EXPECT_TRUE(this->update() == Inertia(1, 1, 0));
   Vector<double> b{1., 2.};
   this->correction->apply(b);
   EXPECT_NEAR(b[0], -1. / 3., 1e-14);
   EXPECT_DOUBLE_EQ(b[1], 2.);
   this->hessian->coefficients = {1.};
   EXPECT_TRUE(this->update() == Inertia(2, 0, 0));
}

TEST_F(LowRankCorrectionTest, PositiveRankOneTermCanRemoveNegativeBaseCurvature) {
   this->solver.diagonal = {1., -1.};
   this->hessian->columns = {{0., 2.}};
   this->hessian->coefficients = {1.};
   EXPECT_TRUE(this->update() == Inertia(2, 0, 0));
}

TEST_F(LowRankCorrectionTest, SingularCorrectionIsReportedToRegularization) {
   this->hessian->columns = {{1., 0.}};
   this->hessian->coefficients = {-1.};
   EXPECT_TRUE(this->update() == Inertia(1, 0, 1));
   Vector<double> b{1., 2.};
   EXPECT_THROW(this->correction->apply(b), std::runtime_error);
}

TEST_F(LowRankCorrectionTest, TwoByTwoPivotHasCorrectInertiaAndSolve) {
   this->hessian->columns = {{1., 0.}, {1., 1.}};
   this->hessian->coefficients = {-1., -0.5};
   // T = [[0, 1], [1, 0]], K = [[-0.5, -0.5], [-0.5, 0.5]].
   EXPECT_TRUE(this->update() == Inertia(1, 1, 0));
   Vector<double> b{1., 2.};
   this->correction->apply(b);
   EXPECT_NEAR(-0.5 * b[0] - 0.5 * b[1], 1., 1e-14);
   EXPECT_NEAR(-0.5 * b[0] + 0.5 * b[1], 2., 1e-14);
}

TEST_F(LowRankCorrectionTest, InvalidCoefficientsAreRejected) {
   this->hessian->columns = {{1., 0.}};
   for (double value: {0., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
      this->hessian->coefficients = {value};
      EXPECT_THROW(this->update(), std::runtime_error);
   }
}

TEST_F(LowRankCorrectionTest, SingularFullMatrixTriggersRegularization) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "requires a direct linear solver";
   }
   this->hessian->columns = {{1., 0.}};
   this->hessian->coefficients = {-1.};
   PrimalDualInertiaCorrection regularization(this->options);
   Subproblem subproblem(this->problem, *this->hessian, regularization);
   Statistics statistics;
   regularization.initialize_statistics(statistics);
   EQPSolver<WoodburyCorrection> solver(this->options, WoodburyCorrection(*this->hessian));
   solver.initialize_memory(subproblem);
   solver.initialize_statistics(statistics);
   Iterate iterate(2, 0);
   this->model.initial_primal_point(iterate.primals);
   Evaluations evaluations(this->model);
   const Direction& direction = solver.solve(statistics, subproblem, iterate, Inf, Vector<double>(2), evaluations, {});
   ASSERT_EQ(direction.status, SubproblemStatus::OPTIMAL);
   const double shift = regularization.get_primal_regularization_factor();
   ASSERT_GT(shift, 0.);
   EXPECT_NEAR(shift * direction.primals[0], -1., 1e-10);
   EXPECT_NEAR((1. + shift) * direction.primals[1], -2., 1e-12);
}
