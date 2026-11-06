// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <climits>
#include <stdexcept>
#include <string>
#include "KrylovSolver.hpp"

namespace uno {
   namespace {
      int to_int(size_t value) {
         if (INT_MAX < value) {
            throw std::overflow_error("libkrylov: dimension too large");
         }
         return static_cast<int>(value);
      }

      // run an operator from a libkrylov callback: an exception must not unwind through the Julia frames, it is
      // stored and rethrown after krylov_solve returns
      template <typename Operators>
      void apply(const KrylovOperator* op, const void* x, void* y, Operators& operators) {
         if (operators.exception) {
            return;
         }
         try {
            (*op)(static_cast<const double*>(x), static_cast<double*>(y));
         }
         catch (...) {
            operators.exception = std::current_exception();
         }
      }
   } // namespace

   KrylovSolver::KrylovSolver(KrylovMethod method, size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns) {
      const KrylovLibrary& library = krylov_library();
      if (library.abi == 0) {
         throw std::runtime_error("libkrylov is not loaded (check the option libkrylov_path or the UNO_KRYLOV_LIBRARY "
            "environment variable)");
      }
      const KrylovWorkspaceOptions workspace_options = library.default_workspace_options();
      const int status = library.workspace_create(static_cast<int>(method), to_int(number_rows), to_int(number_columns),
         KRYLOV_FLOAT64, KRYLOV_CPU, &workspace_options, &this->workspace);
      if (status != 0 || this->workspace == nullptr) {
         throw std::runtime_error("libkrylov: krylov_workspace_create failed with status " + std::to_string(status));
      }
   }

   KrylovSolver::~KrylovSolver() {
      if (this->workspace != nullptr) {
         krylov_library().workspace_free(this->workspace);
      }
   }

   void KrylovSolver::apply_A(const void* x, void* y, void* userdata) {
      auto& operators = *static_cast<Operators*>(userdata);
      apply(operators.A, x, y, operators);
   }

   void KrylovSolver::apply_At(const void* x, void* y, void* userdata) {
      auto& operators = *static_cast<Operators*>(userdata);
      apply(operators.At, x, y, operators);
   }

   void KrylovSolver::apply_M(const void* x, void* y, void* userdata) {
      auto& operators = *static_cast<Operators*>(userdata);
      apply(operators.M, x, y, operators);
   }

   KrylovResult KrylovSolver::solve(const KrylovOperator& A, const double* b, double* x, const KrylovParameters& parameters,
         const KrylovOperator* At, const KrylovOperator* M) {
      const KrylovLibrary& library = krylov_library();
      Operators operators{&A, At, M, nullptr};
      const KrylovMatvec matvec_At = (At != nullptr) ? &KrylovSolver::apply_At : nullptr;
      const KrylovMatvec matvec_M = (M != nullptr) ? &KrylovSolver::apply_M : nullptr;

      int status;
      if (library.abi == 2) {
         KrylovOptionsABI2 options = library.default_options_abi2();
         options.atol = parameters.atol;
         options.rtol = parameters.rtol;
         options.itmax = parameters.itmax;
         options.verbose = parameters.verbose;
         options.lambda = parameters.lambda;
         status = library.solve_abi2(this->workspace, &KrylovSolver::apply_A, matvec_At, matvec_M, nullptr, b, nullptr,
            &operators, &options);
      }
      else {
         KrylovOptionsABI1 options = library.default_options_abi1();
         options.atol = parameters.atol;
         options.rtol = parameters.rtol;
         options.itmax = parameters.itmax;
         options.verbose = parameters.verbose;
         options.lambda = parameters.lambda;
         status = library.solve_abi1(this->workspace, &KrylovSolver::apply_A, matvec_At, matvec_M, b, nullptr,
            &operators, &options);
      }
      if (operators.exception) {
         std::rethrow_exception(operators.exception);
      }
      if (status != 0) {
         throw std::runtime_error("libkrylov: krylov_solve failed with status " + std::to_string(status));
      }
      if (library.get_x(this->workspace, x, to_int(this->number_columns)) != 0) {
         throw std::runtime_error("libkrylov: krylov_get_x failed");
      }
      const int iterations = library.niter(this->workspace);
      return {library.is_solved(this->workspace) == 1, static_cast<size_t>(iterations < 0 ? 0 : iterations),
         library.elapsed_time(this->workspace)};
   }

   bool KrylovSolver::warm_start(const double* x0) {
      return krylov_library().warm_start(this->workspace, x0, to_int(this->number_columns)) == 0;
   }
} // namespace
