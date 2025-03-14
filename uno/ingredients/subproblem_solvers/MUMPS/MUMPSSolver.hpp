// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MUMPSSOLVER_H
#define UNO_MUMPSSOLVER_H

#include <vector>
#include "dmumps_c.h"
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class MUMPSSolver : public DirectSymmetricIndefiniteLinearSolver<size_t, double> {
   public:
      explicit MUMPSSolver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros);
      ~MUMPSSolver() override;

      void do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) override;
      void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
      void solve_EQP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, Vector<double>& result,
            WarmstartInformation& warmstart_information) override;

      [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] size_t number_zero_eigenvalues() const;
      // [[nodiscard]] bool matrix_is_positive_definite() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

   protected:
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */
      SymmetricMatrix<size_t, double> hessian;

      const size_t dimension;
      const size_t number_nonzeros;

      DMUMPS_STRUC_C mumps_structure{};

      // matrix sparsity
      std::vector<int> row_indices{};
      std::vector<int> column_indices{};
      SymmetricMatrix<size_t, double> augmented_matrix;
      Vector<double> rhs;

      static const int JOB_INIT = -1;
      static const int JOB_END = -2;
      static const int JOB_ANALYSIS = 1;
      static const int JOB_FACTORIZATION = 2;
      static const int JOB_SOLVE = 3;

      static const int GENERAL_SYMMETRIC = 2;

      const size_t fortran_shift{1};
      void set_up_subproblem(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, WarmstartInformation& warmstart_information);
      void save_sparsity_to_local_format(const SymmetricMatrix<size_t, double>& matrix);
      void solve_indefinite_linear_system(Vector<double>& result);
   };
} // namespace

#endif // UNO_MUMPSSOLVER_H
