// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SSIDSSOLVER_H
#define UNO_SSIDSSOLVER_H

#include "spral_ssids.h"
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   struct Workspace {
      void* akeep{nullptr};
      void* fkeep{nullptr};
      spral_ssids_options options{};
      spral_ssids_inform inform{};
      int n;
      int nnz;
   };

   class SSIDSSolver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      SSIDSSolver();
      ~SSIDSSolver() override = default;

      void initialize_memory() override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* result) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system();

   protected:
      Workspace workspace{};
      COOLinearSystem linear_system{Indexing::Fortran_indexing};

      bool analysis_performed{false};
      bool factorization_performed{false};
   };
} // namespace

#endif // UNO_SSIDSSOLVER_H