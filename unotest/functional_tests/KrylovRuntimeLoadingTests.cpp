// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Runtime loading of libkrylov. The library is resolved through the UNO_KRYLOV_LIBRARY environment variable
// (e.g. dependencies/libkrylov/lib/libkrylov.so); the tests are skipped if it is not set.

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/Krylov/KrylovLoader.hpp"
#include "ingredients/subproblem_solvers/Krylov/KrylovSolver.hpp"

using namespace uno;

namespace {
   bool libkrylov_requested() {
      const char* path = std::getenv("UNO_KRYLOV_LIBRARY");
      return path != nullptr && path[0] != '\0';
   }
}

TEST(KrylovRuntimeLoading, IsFunctional) {
   if (!libkrylov_requested()) {
      GTEST_SKIP() << "UNO_KRYLOV_LIBRARY is not set";
   }
   ASSERT_TRUE(is_libkrylov_functional()) << krylov_library_error();
   ASSERT_NE(krylov_library().abi, 0);
   ASSERT_FALSE(loaded_krylov_library().empty());
}

TEST(KrylovRuntimeLoading, MINRESOnIndefiniteSystem) {
   if (!libkrylov_requested() || !load_krylov_library()) {
      GTEST_SKIP() << "libkrylov is not available";
   }
   // symmetric indefinite (KKT-like) system [[2, 1], [1, -1]] x = [3, 0] -> x = (1, 1)
   const KrylovOperator A = [](const double* x, double* y) {
      y[0] = 2. * x[0] + x[1];
      y[1] = x[0] - x[1];
   };
   const double b[2] = {3., 0.};
   double x[2] = {0., 0.};
   KrylovSolver solver(KrylovMethod::MINRES, 2, 2);
   KrylovParameters parameters;
   parameters.atol = parameters.rtol = 1e-12;
   const KrylovResult result = solver.solve(A, b, x, parameters);
   ASSERT_TRUE(result.solved);
   ASSERT_NEAR(x[0], 1., 1e-8);
   ASSERT_NEAR(x[1], 1., 1e-8);
}

TEST(KrylovRuntimeLoading, LSMRWithAdjoint) {
   if (!libkrylov_requested() || !load_krylov_library()) {
      GTEST_SKIP() << "libkrylov is not available";
   }
   // overdetermined least-squares problem min ||A x - b||, A = [[1, 0], [0, 1], [1, 1]], b = [1, 1, 2] -> x = (1, 1)
   const KrylovOperator A = [](const double* x, double* y) {
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[0] + x[1];
   };
   const KrylovOperator At = [](const double* y, double* x) {
      x[0] = y[0] + y[2];
      x[1] = y[1] + y[2];
   };
   const double b[3] = {1., 1., 2.};
   double x[2] = {0., 0.};
   KrylovSolver solver(KrylovMethod::LSMR, 3, 2);
   KrylovParameters parameters;
   parameters.atol = parameters.rtol = 1e-12;
   const KrylovResult result = solver.solve(A, b, x, parameters, &At);
   ASSERT_TRUE(result.solved);
   ASSERT_NEAR(x[0], 1., 1e-8);
   ASSERT_NEAR(x[1], 1., 1e-8);
}

TEST(KrylovRuntimeLoading, CallbackExceptionIsRethrown) {
   if (!libkrylov_requested() || !load_krylov_library()) {
      GTEST_SKIP() << "libkrylov is not available";
   }
   const KrylovOperator A = [](const double* /*x*/, double* /*y*/) {
      throw std::runtime_error("evaluation error");
   };
   const double b[2] = {1., 1.};
   double x[2] = {0., 0.};
   KrylovSolver solver(KrylovMethod::CG, 2, 2);
   ASSERT_THROW((void) solver.solve(A, b, x, KrylovParameters{}), std::runtime_error);
}
