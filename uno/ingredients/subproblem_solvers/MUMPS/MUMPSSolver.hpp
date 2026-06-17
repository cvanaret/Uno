// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MUMPSSOLVER_H
#define UNO_MUMPSSOLVER_H

#include <dmumps_c.h>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // settings taken from MadNLP
   struct MUMPSSettings {
      static constexpr int mem_percent = 35;
      static constexpr int mem_percent_increase = 2;
      static constexpr int permuting_scaling = 0;
      static constexpr int pivot_order = 0; // AMD
      static constexpr double pivtol = 1e-6; // relative threshold for numerical pivoting
      static constexpr double pivtolmax = 1e-1;
      static constexpr int scaling = 1; // diagonal scaling computed during the numerical factorization phase
      static constexpr size_t max_number_factorization_failures = 10;
   };

   class MUMPSSolver : public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      MUMPSSolver();
      ~MUMPSSolver() override;

      void initialize_memory() override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* result) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] size_t number_zero_eigenvalues() const;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system();

   protected:
      DMUMPS_STRUC_C workspace{};
      COOLinearSystem linear_system{Indexing::Fortran_indexing};

      static const int JOB_INIT = -1;
      static const int JOB_END = -2;
      static const int JOB_ANALYSIS = 1;
      static const int JOB_FACTORIZATION = 2;
      static const int JOB_SOLVE = 3;

      static const int GENERAL_SYMMETRIC = 2;

      bool analysis_performed{false};
      bool factorization_performed{false};
      size_t number_factorization_failures{0};
      const size_t max_number_factorization_failures{MUMPSSettings::max_number_factorization_failures};

      [[nodiscard]] int& ICNTL(size_t index);
      [[nodiscard]] double& CNTL(size_t index);
      [[nodiscard]] int INFO(size_t index) const;
      [[nodiscard]] int INFOG(size_t index) const;
   };
} // namespace

#endif // UNO_MUMPSSOLVER_H