// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Models that provide a Lagrangian Hessian operator but no explicit Hessian matrix

#include <string>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // hs015 with a Hessian operator only (the dense lower triangle is used to compute the products)
   class HessianOperatorOnlyModel: public DenseTestModel {
   public:
      HessianOperatorOnlyModel(): DenseTestModel("hs015", 2, 2, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return 100.*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1. - x[0])*(1. - x[0]); },
            [](const Vector<double>& x, Vector<double>& g) {
               g[0] = 400.*x[0]*x[0]*x[0] - 400.*x[0]*x[1] + 2.*x[0] - 2.;
               g[1] = 200.*(x[1] - x[0]*x[0]);
            },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[1]; c[1] = x[0] + x[1]*x[1]; },
            [](const Vector<double>& x, double* J) { J[0] = x[1]; J[1] = x[0]; J[2] = 1.; J[3] = 2.*x[1]; },
            [](const Vector<double>& x, double rho, const Vector<double>& y, double* H) {
               H[0] = rho*(1200.*x[0]*x[0] - 400.*x[1] + 2.);
               H[1] = -400.*rho*x[0] - y[0];
               H[2] = 200.*rho - 2.*y[1];
            },
            {-Inf, -Inf}, {0.5, Inf}, {1., 0.}, {Inf, Inf}, {-2., 1.}) {
         this->enable_hessian_operator();
      }

      [[nodiscard]] bool has_hessian_matrix() const override { return false; }
      void evaluate_lagrangian_hessian(const Vector<double>&, double, const Vector<double>&, View<double>) const override {
         ADD_FAILURE() << "the explicit Hessian matrix was evaluated";
      }
   };

   Options preset_options(const std::string& preset) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      Logger::level = SILENT;
      return options;
   }

   bool has_lbfgs_override(const Result& result) {
      for (const OptionOverride& option_override: result.option_overrides) {
         if (option_override.option_name == "hessian_model" && option_override.new_value == "LBFGS") {
            return true;
         }
      }
      return false;
   }

   void check_hs015_solution(const Result& result) {
      EXPECT_NEAR(result.primal_solution[0], 0.5, 1e-5);
      EXPECT_NEAR(result.primal_solution[1], 2., 1e-5);
      EXPECT_NEAR(result.solution_objective, 306.5, 1e-3);
   }
} // namespace

// Bug: with a Hessian operator and no Hessian matrix, the IPOPT preset failed in EQPSolver::initialize_memory (the direct
// linear solver needs the matrix) with ALGORITHMIC_ERROR and 0 iterations. Uno now detects it when the ingredients are
// allocated and falls back to L-BFGS, as when no Hessian is provided
TEST(HessianOperatorOnly, InteriorPointFallsBackToLBFGS) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const HessianOperatorOnlyModel model;
   Uno uno;
   const Result result = uno.solve(model, preset_options("ipopt"));
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
   EXPECT_GT(result.number_iterations, 0);
   EXPECT_TRUE(has_lbfgs_override(result));
   EXPECT_NE(uno.get_method_description().find("L-BFGS"), std::string::npos) << uno.get_method_description();
   EXPECT_EQ(model.number_model_hessian_vector_products(), 0);
   check_hs015_solution(result);
}

// the operator is still used when the subproblem solver can use it (BQPD)
TEST(HessianOperatorOnly, FilterSQPWithBQPDUsesTheOperator) {
#ifndef HAS_BQPD
   GTEST_SKIP() << "the test requires BQPD";
#else
   const HessianOperatorOnlyModel model;
   Options options = preset_options("filtersqp");
   options.set_string("QP_solver", "BQPD");
   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_FALSE(has_lbfgs_override(result));
   EXPECT_NE(uno.get_method_description().find("exact Hessian"), std::string::npos) << uno.get_method_description();
   EXPECT_GT(model.number_model_hessian_vector_products(), 0);
   check_hs015_solution(result);
#endif
}
