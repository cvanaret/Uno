// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_KRYLOVSOLVER_H
#define UNO_KRYLOVSOLVER_H

#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include "KrylovLoader.hpp"

namespace uno {
   // y = A x (or A^T x, or y = M^{-1} x for a preconditioner)
   using KrylovOperator = std::function<void(const double* x, double* y)>;

   // solve-time parameters. NaN/0 sentinels: libkrylov defaults
   struct KrylovParameters {
      double atol{std::numeric_limits<double>::quiet_NaN()};
      double rtol{std::numeric_limits<double>::quiet_NaN()};
      int itmax{0};
      int verbose{0};
      double lambda{0.}; // regularization/shift
   };

   struct KrylovResult {
      bool solved;
      size_t iterations;
      double time;
   };

   // RAII wrapper around a libkrylov workspace (double precision, CPU) for an m x n operator.
   // libkrylov must have been loaded with load_krylov_library()
   class KrylovSolver {
   public:
      KrylovSolver(KrylovMethod method, size_t number_rows, size_t number_columns);
      ~KrylovSolver();
      KrylovSolver(const KrylovSolver&) = delete;
      KrylovSolver& operator=(const KrylovSolver&) = delete;

      // solve with the operator A, right-hand side b (size m) and solution x (size n). The adjoint A^T (required by
      // e.g. LSQR, LSMR, CGLS, BiLQ, QMR) and the preconditioner M are optional
      KrylovResult solve(const KrylovOperator& A, const double* b, double* x, const KrylovParameters& parameters,
         const KrylovOperator* At = nullptr, const KrylovOperator* M = nullptr);
      // initial guess (size n) for the next solve. Returns false if the method does not support warm starting
      bool warm_start(const double* x0);

   protected:
      const size_t number_rows;
      const size_t number_columns;
      void* workspace{nullptr};

      // passed as userdata to the trampolines
      struct Operators {
         const KrylovOperator* A;
         const KrylovOperator* At;
         const KrylovOperator* M;
         std::exception_ptr exception;
      };
      static void apply_A(const void* x, void* y, void* userdata);
      static void apply_At(const void* x, void* y, void* userdata);
      static void apply_M(const void* x, void* y, void* userdata);
   };
} // namespace

#endif // UNO_KRYLOVSOLVER_H
