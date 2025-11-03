// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_USERMODEL_H
#define UNO_USERMODEL_H

#include <cstring>
#include <optional>
#include <vector>
#include "C/Uno_C_API.h"
#include "optimization/ProblemType.hpp"

namespace uno {
   inline ProblemType problem_type_from_string(const char* problem_type) {
      if (strcmp(problem_type, UNO_PROBLEM_LINEAR) == 0) {
         return ProblemType::LINEAR;
      }
      else if (strcmp(problem_type, UNO_PROBLEM_QUADRATIC) == 0) {
         return ProblemType::QUADRATIC;
      }
      else {
         return ProblemType::NONLINEAR;
      }
   }

   // UserModel contains the description of the model provided by the user
   template <typename Objective, typename ObjectiveGradient, typename Constraints, typename Jacobian,
      typename JacobianOperator, typename JacobianTransposedOperator, typename Hessian, typename HessianOperator,
      typename DoubleVector, typename UserDataType>
   class UserModel {
   public:
      // problem_type is "LP" for linear, "QP" for quadratic, "NLP" for nonlinear
      UserModel(const char* problem_type, uno_int number_variables, uno_int base_indexing):
            problem_type(problem_type_from_string(problem_type)),
            base_indexing(base_indexing),
            number_variables(number_variables) {
      }

      ~UserModel() = default;

      const ProblemType problem_type;
      const uno_int base_indexing; // 0 for C-style indexing, 1 for Fortran-style indexing

      // variables
      const uno_int number_variables;
      DoubleVector variables_lower_bounds{};
      DoubleVector variables_upper_bounds{};

      // objective
      Objective objective_function{nullptr};
      ObjectiveGradient objective_gradient{nullptr};

      // constraints
      uno_int number_constraints{0};
      Constraints constraint_functions{nullptr};
      DoubleVector constraints_lower_bounds{};
      DoubleVector constraints_upper_bounds{};
      uno_int number_jacobian_nonzeros{0};
      std::vector<uno_int> jacobian_row_indices{};
      std::vector<uno_int> jacobian_column_indices{};
      Jacobian constraint_jacobian{nullptr};
      JacobianOperator jacobian_operator{nullptr};
      JacobianTransposedOperator jacobian_transposed_operator{nullptr};

      // Hessian
      std::optional<uno_int> number_hessian_nonzeros{};
      // lower ('L') or upper ('U')
      char hessian_triangular_part{}; // default is empty
      std::vector<uno_int> hessian_row_indices{};
      std::vector<uno_int> hessian_column_indices{};
      Hessian lagrangian_hessian{nullptr};
      HessianOperator lagrangian_hessian_operator{nullptr};
      double lagrangian_sign_convention{UNO_MULTIPLIER_NEGATIVE};

      // User data
      UserDataType user_data{};

      // Optimization sense
      uno_int optimization_sense{UNO_MINIMIZE};

      // initial iterate
      DoubleVector initial_primal_iterate{};
      DoubleVector initial_dual_iterate{};
   };
} // namespace

#endif // UNO_USERMODEL_H