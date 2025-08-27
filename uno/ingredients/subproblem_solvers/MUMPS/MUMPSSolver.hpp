// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MUMPSSOLVER_H
#define UNO_MUMPSSOLVER_H

#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "dmumps_c.h"
#include "../COOEvaluationSpace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class MUMPSSolver : public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      MUMPSSolver();
      ~MUMPSSolver() override;

      void initialize_hessian(const Subproblem& subproblem) override;
      void initialize_augmented_system(const Subproblem& subproblem) override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(const double* matrix_values) override;
      void solve_indefinite_system(const Vector<double>& matrix_values, const Vector<double>& rhs, Vector<double>& result) override;
      void solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] size_t number_zero_eigenvalues() const;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override;

   protected:
      DMUMPS_STRUC_C workspace{};
      COOEvaluationSpace evaluation_space{};

      static const int JOB_INIT = -1;
      static const int JOB_END = -2;
      static const int JOB_ANALYSIS = 1;
      static const int JOB_FACTORIZATION = 2;
      static const int JOB_SOLVE = 3;

      static const int GENERAL_SYMMETRIC = 2;

      bool analysis_performed{false};
      bool factorization_performed{false};
   };
} // namespace

#endif // UNO_MUMPSSOLVER_H