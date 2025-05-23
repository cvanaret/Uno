// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MUMPSSOLVER_H
#define UNO_MUMPSSOLVER_H

#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "dmumps_c.h"

namespace uno {
   class MUMPSSolver : public DirectSymmetricIndefiniteLinearSolver<size_t, double> {
   public:
      MUMPSSolver();
      ~MUMPSSolver() override;

      void initialize_memory(size_t dimension, size_t number_nonzeros) override;

      void do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) override;
      void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
      void solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] size_t number_zero_eigenvalues() const;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

   protected:
      DMUMPS_STRUC_C mumps_structure{};
      size_t dimension{};

      // matrix sparsity
      std::vector<int> row_indices{};
      std::vector<int> column_indices{};

      static const int JOB_INIT = -1;
      static const int JOB_END = -2;
      static const int JOB_ANALYSIS = 1;
      static const int JOB_FACTORIZATION = 2;
      static const int JOB_SOLVE = 3;

      static const int GENERAL_SYMMETRIC = 2;

      const size_t fortran_shift{1};
      void save_sparsity_to_local_format(const SymmetricMatrix<size_t, double>& matrix);
   };
} // namespace

#endif // UNO_MUMPSSOLVER_H
