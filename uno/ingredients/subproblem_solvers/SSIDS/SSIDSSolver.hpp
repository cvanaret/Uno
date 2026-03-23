// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SSIDSSOLVER_H
#define UNO_SSIDSSOLVER_H

#include <spral_ssids.h>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOWorkspace.hpp"

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

      void initialize_hessian(const Subproblem& subproblem) override;
      void initialize_augmented_system(const Subproblem& subproblem) override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(const double* matrix_values, bool is_matrix_positive_definite) override;
      void solve_indefinite_system(const double* matrix_values, const double* rhs, double* result) override;
      void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      Workspace workspace{};
      COOWorkspace coo_workspace{};

      bool analysis_performed{false};
      bool factorization_performed{false};
   };
} // namespace

#endif // UNO_SSIDSSOLVER_H