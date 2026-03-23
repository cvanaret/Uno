// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOWORKSPACE_H
#define UNO_COOWORKSPACE_H

#include <cstddef>
#include <vector>
#include "LinearSolverSparseRepresentation.hpp"
#include "SolverWorkspace.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Evaluations;
   class Statistics;
   class WarmstartInformation;

   class COOLinearSolverSparseRepresentation: public LinearSolverSparseRepresentation {
   public:
      COOLinearSolverSparseRepresentation() = default;
      ~COOLinearSolverSparseRepresentation() override = default;

      void initialize_hessian(const Subproblem& subproblem) override;
      void initialize_augmented_system(const Subproblem& subproblem) override;

      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const override;

      Vector<double> objective_gradient{}; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints{}; /*!< Constraint values (size \f$m)\f$ */

      // Jacobian
      size_t number_jacobian_nonzeros{};

      // symmetric matrix (Hessian or augmented system)
      size_t number_hessian_nonzeros{};
      std::vector<uno_int> matrix_row_indices{};
      std::vector<uno_int> matrix_column_indices{};
      bool analysis_performed{false};
   };
} // namespace

#endif // UNO_COOWORKSPACE_H