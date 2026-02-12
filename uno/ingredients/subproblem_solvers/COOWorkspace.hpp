// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOWORKSPACE_H
#define UNO_COOWORKSPACE_H

#include <cstddef>
#include <vector>
#include "linear_algebra/Vector.hpp"
#include "SolverWorkspace.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   class COOWorkspace: public SolverWorkspace {
   public:
      COOWorkspace() = default;
      ~COOWorkspace() override = default;

      void initialize_hessian(const Subproblem& subproblem);
      void initialize_augmented_system(const Subproblem& subproblem);

      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void set_up_linear_system(Statistics& statistics, const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         Evaluations& evaluations, const WarmstartInformation& warmstart_information);

      Vector<double> objective_gradient{}; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints{}; /*!< Constraint values (size \f$m)\f$ */

      // Jacobian
      size_t number_jacobian_nonzeros{};
      std::vector<uno_int> jacobian_row_indices{};
      std::vector<uno_int> jacobian_column_indices{};

      // symmetric matrix (Hessian or augmented system)
      size_t number_hessian_nonzeros{};
      size_t number_matrix_nonzeros{};
      std::vector<uno_int> matrix_row_indices{};
      std::vector<uno_int> matrix_column_indices{};
      Vector<double> matrix_values;
      Vector<double> rhs{};
      Vector<double> solution{};
      bool analysis_performed{false};
   };
} // namespace

#endif // UNO_COOWORKSPACE_H