// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_USERMODEL_H
#define UNO_USERMODEL_H

#include <vector>
#include "C/Uno_C_API.h"

namespace uno {
   // UserModel contains the description of the model provided by the user
   template <typename Objective, typename ObjectiveGradient, typename Constraints, typename Jacobian,
      typename JacobianOperator, typename JacobianTransposedOperator, typename Hessian, typename HessianOperator,
      typename DoubleVector, typename UserDataType>
   class UserModel {
   public:
      UserModel(char problem_type, int32_t number_variables, DoubleVector variables_lower_bounds,
         DoubleVector variables_upper_bounds, int32_t base_indexing):
            problem_type(problem_type),
            base_indexing(base_indexing),
            number_variables(number_variables),
            variables_lower_bounds(variables_lower_bounds),
            variables_upper_bounds(variables_upper_bounds) {
      }

      ~UserModel() = default;

      const char problem_type; // 'L' for linear, 'Q' for quadratic, 'N' for nonlinear
      const int32_t base_indexing; // 0 for C-style indexing, 1 for Fortran-style indexing

      // variables
      const int32_t number_variables;
      DoubleVector variables_lower_bounds{nullptr};
      DoubleVector variables_upper_bounds{nullptr};

      // objective
      Objective objective_function{nullptr};
      ObjectiveGradient objective_gradient{nullptr};

      // constraints
      int32_t number_constraints{0};
      Constraints constraint_functions{nullptr};
      DoubleVector constraints_lower_bounds{};
      DoubleVector constraints_upper_bounds{};
      int32_t number_jacobian_nonzeros{0};
      std::vector<int32_t> jacobian_row_indices{};
      std::vector<int32_t> jacobian_column_indices{};
      Jacobian constraint_jacobian{nullptr};
      JacobianOperator jacobian_operator{nullptr};
      JacobianTransposedOperator jacobian_transposed_operator{nullptr};

      // Hessian
      int32_t number_hessian_nonzeros{0};
      // lower ('L') or upper ('U')
      char hessian_triangular_part{}; // default is empty
      std::vector<int32_t> hessian_row_indices{};
      std::vector<int32_t> hessian_column_indices{};
      Hessian lagrangian_hessian{nullptr};
      HessianOperator lagrangian_hessian_operator{nullptr};
      double lagrangian_sign_convention{UNO_MULTIPLIER_NEGATIVE};

      // User data
      UserDataType user_data{};

      // Optimization sense
      int32_t optimization_sense{UNO_MINIMIZE};

      // initial iterate
      DoubleVector initial_primal_iterate{};
      DoubleVector initial_dual_iterate{};
   };
} // namespace

#endif // UNO_USERMODEL_H