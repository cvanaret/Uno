// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Failures of the Hessian-vector product called by BQPD through the Fortran gdotx callback

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/hessian_models/ExactHessian.hpp"
#include "ingredients/inertia_correction_strategies/NoInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/BQPD/BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/BQPD/BQPDSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/Result.hpp"
#include "optimization/UserTermination.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "fortran_interface.h"

using namespace uno;

// BQPD's Hessian-vector product callback, defined in BQPDSolver.cpp
extern "C" void FC_GLOBAL(gdotx, GDOTX)(int* dimension, const double vector[], const double ws[], const int lws[], double result[]);

namespace {
   enum class HessianFailure { NONE, THROW, NAN_PRODUCT, USER_TERMINATION, RUNTIME_ERROR };

   // min x0^2 + x1^2  s.t.  x0 + x1 <= 1.5, -10 <= x <= 10, with a Hessian operator that fails at its failing_call-th
   // call (1-based), either by throwing a HessianEvaluationError or by returning NaNs
   class FailingHessianOperatorModel: public DenseTestModel {
   public:
      FailingHessianOperatorModel(): DenseTestModel("failing_hessian_operator", 2, 1, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1]; },
            [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
            [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
            [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
            {-10., -10.}, {10., 10.}, {-Inf}, {1.5}, {0.3, 0.4}) {
         this->enable_hessian_operator();
      }

      // the operator only (no explicit Hessian matrix): BQPD uses the gdotx callback
      [[nodiscard]] bool has_hessian_matrix() const override { return false; }

      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
            const Vector<double>& multipliers, View<double> result) const override {
         ++this->number_operator_calls;
         if (this->number_operator_calls == this->failing_call) {
            if (this->failure == HessianFailure::THROW) {
               throw HessianEvaluationError();
            }
            if (this->failure == HessianFailure::USER_TERMINATION) {
               throw UserTermination();
            }
            if (this->failure == HessianFailure::RUNTIME_ERROR) {
               throw std::runtime_error("operator failure");
            }
            if (this->failure == HessianFailure::NAN_PRODUCT) {
               result.fill(std::numeric_limits<double>::quiet_NaN());
               return;
            }
         }
         DenseTestModel::compute_hessian_vector_product(x, vector, objective_multiplier, multipliers, result);
      }

      HessianFailure failure{HessianFailure::NONE};
      size_t failing_call{1};
      mutable size_t number_operator_calls{0};
   };

   // solves the QP subproblem at (0.3, 0.4) with BQPD: min g^T d + d^T d (g = 2x)  s.t.  d0 + d1 <= 0.8
   class BQPDSubproblemFixture {
   public:
      explicit BQPDSubproblemFixture(const FailingHessianOperatorModel& model):
            problem(model, 100.),
            hessian_model(model),
            subproblem(this->problem, this->hessian_model, this->inertia_correction_strategy),
            solver(BQPDSubproblemFixture::options()),
            iterate(this->subproblem.number_variables, this->subproblem.number_constraints),
            evaluations(model),
            direction(this->subproblem.number_variables, this->subproblem.number_constraints),
            initial_point(this->subproblem.number_variables) {
         Logger::level = SILENT;
         this->solver.initialize_memory(this->subproblem);
         this->iterate.primals[0] = 0.3;
         this->iterate.primals[1] = 0.4;
         this->solver.get_quadratic_program().fill(this->statistics, this->subproblem, this->iterate,
            /* trust_region_radius = */ Inf, this->evaluations, this->warmstart_information);
      }

      void solve() {
         this->solver.solve(this->statistics, this->initial_point, this->direction, this->warmstart_information);
      }

      static Options options() {
         Options options;
         DefaultOptions::load(options);
         options.set_bool("print_subproblem", false);
         return options;
      }

      OptimizationProblem problem;
      ExactHessian hessian_model;
      NoInertiaCorrection inertia_correction_strategy{};
      Subproblem subproblem;
      BQPDSolver solver;
      Iterate iterate;
      Evaluations evaluations;
      Direction direction;
      Vector<double> initial_point;
      WarmstartInformation warmstart_information{};
      Statistics statistics{};
   };

   Options filtersqp_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "filtersqp");
      options.set_string("QP_solver", "BQPD");
      Logger::level = SILENT;
      return options;
   }
} // namespace

// Bug: a Hessian operator that threw inside BQPD's gdotx callback unwound through the Fortran frames of BQPD
// (undefined behavior). The exception is now stored, BQPD receives zero products without calling the operator again,
// and the exception is rethrown once BQPD has returned. The next solve starts afresh.
TEST(BQPDHessianErrors, OperatorExceptionIsRethrownAfterBQPDReturns) {
   FailingHessianOperatorModel model;
   model.failure = HessianFailure::THROW;
   BQPDSubproblemFixture fixture(model);

   EXPECT_THROW(fixture.solve(), HessianEvaluationError);
   EXPECT_EQ(model.number_operator_calls, 1) << "the operator was called again after its failure";

   // the operator no longer fails: the same subproblem is solved correctly (d = -x)
   const size_t number_calls_before = model.number_operator_calls;
   fixture.solve();
   EXPECT_GT(model.number_operator_calls, number_calls_before);
   ASSERT_EQ(fixture.direction.status, SubproblemStatus::OPTIMAL);
   EXPECT_NEAR(fixture.direction.primals[0], -0.3, 1e-8);
   EXPECT_NEAR(fixture.direction.primals[1], -0.4, 1e-8);
}

// the gdotx callback (called by the Fortran code of BQPD) never lets an exception escape: it returns zero products and
// defers the exception
TEST(BQPDHessianErrors, GdotxCallbackDoesNotThrow) {
   FailingHessianOperatorModel model;
   model.failure = HessianFailure::THROW;
   BQPDSubproblemFixture fixture(model);
   auto& quadratic_program = dynamic_cast<BQPDQuadraticProgram&>(fixture.solver.get_quadratic_program());
   std::vector<int> lws(hidden_pointer_size_in_lws);
   hide_pointer(0, lws.data(), quadratic_program);

   int dimension = 2;
   const double vector[2] = {1., 1.};
   const double ws[1] = {0.};
   double result[2] = {7., 7.};
   EXPECT_NO_THROW(FC_GLOBAL(gdotx, GDOTX)(&dimension, vector, ws, lws.data(), result));
   EXPECT_EQ(result[0], 0.);
   EXPECT_EQ(result[1], 0.);
   // the subsequent calls return zero products without calling the operator
   result[0] = result[1] = 7.;
   EXPECT_NO_THROW(FC_GLOBAL(gdotx, GDOTX)(&dimension, vector, ws, lws.data(), result));
   EXPECT_EQ(result[0], 0.);
   EXPECT_EQ(result[1], 0.);
   EXPECT_EQ(model.number_operator_calls, 1);
   // the exception is rethrown once, then cleared
   EXPECT_THROW(quadratic_program.rethrow_hessian_vector_product_exception(), HessianEvaluationError);
   EXPECT_NO_THROW(quadratic_program.rethrow_hessian_vector_product_exception());
   // the operator is called again
   EXPECT_NO_THROW(FC_GLOBAL(gdotx, GDOTX)(&dimension, vector, ws, lws.data(), result));
   EXPECT_EQ(model.number_operator_calls, 2);
   EXPECT_DOUBLE_EQ(result[0], 2.);
   EXPECT_DOUBLE_EQ(result[1], 2.);
}

// Bug: a Hessian operator that returned NaNs could make BQPD hang. A non-finite product is now reported as a
// HessianEvaluationError once BQPD has returned
TEST(BQPDHessianErrors, NonFiniteProductIsAnEvaluationError) {
   FailingHessianOperatorModel model;
   model.failure = HessianFailure::NAN_PRODUCT;
   BQPDSubproblemFixture fixture(model);

   EXPECT_THROW(fixture.solve(), HessianEvaluationError);
   EXPECT_EQ(model.number_operator_calls, 1) << "the operator was called again after its failure";

   fixture.solve();
   ASSERT_EQ(fixture.direction.status, SubproblemStatus::OPTIMAL);
   EXPECT_NEAR(fixture.direction.primals[0], -0.3, 1e-8);
   EXPECT_NEAR(fixture.direction.primals[1], -0.4, 1e-8);
}

// end to end: a failing Hessian operator is handled by the trust region like any evaluation error (the radius is
// decreased) and Uno converges
TEST(BQPDHessianErrors, FilterSQPRecoversFromFailingOperator) {
   for (const HessianFailure failure: {HessianFailure::THROW, HessianFailure::NAN_PRODUCT}) {
      FailingHessianOperatorModel model;
      model.failure = failure;
      const Options options = filtersqp_options();
      Uno uno;
      const Result result = uno.solve(model, options);
      const std::string label = (failure == HessianFailure::THROW) ? "throw" : "NaN";
      ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << label;
      ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT) << label;
      EXPECT_NEAR(result.primal_solution[0], 0., 1e-6) << label;
      EXPECT_NEAR(result.primal_solution[1], 0., 1e-6) << label;
      EXPECT_GT(model.number_operator_calls, 1) << label;
   }
}

// the deferred failure keeps its kind: a termination request and other errors are rethrown as such once BQPD has
// returned (the failure is not stored as a std::exception_ptr, which breaks when two C++ runtimes are loaded)
TEST(BQPDHessianErrors, DeferredFailureKeepsItsKind) {
   {
      FailingHessianOperatorModel model;
      model.failure = HessianFailure::USER_TERMINATION;
      BQPDSubproblemFixture fixture(model);
      EXPECT_THROW(fixture.solve(), UserTermination);
      EXPECT_EQ(model.number_operator_calls, 1);
   }
   {
      FailingHessianOperatorModel model;
      model.failure = HessianFailure::RUNTIME_ERROR;
      BQPDSubproblemFixture fixture(model);
      try {
         fixture.solve();
         ADD_FAILURE() << "no exception was rethrown";
      }
      catch (const std::runtime_error& exception) {
         EXPECT_EQ(std::string(exception.what()), "operator failure");
      }
      EXPECT_EQ(model.number_operator_calls, 1);
   }
}

// end to end: a termination request from the Hessian operator stops Uno with USER_TERMINATION
TEST(BQPDHessianErrors, UserTerminationFromOperatorStopsUno) {
   FailingHessianOperatorModel model;
   model.failure = HessianFailure::USER_TERMINATION;
   Uno uno;
   const Result result = uno.solve(model, filtersqp_options());
   EXPECT_EQ(result.optimization_status, OptimizationStatus::USER_TERMINATION);
   EXPECT_EQ(model.number_operator_calls, 1);
}
