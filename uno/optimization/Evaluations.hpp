// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONS_H
#define UNO_EVALUATIONS_H

#include <vector>
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "tools/Infinity.hpp"

struct Evaluations {
   double objective{INF<double>}; /*!< Objective value */
   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
   RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

   Evaluations(size_t max_number_variables, size_t max_number_constraints):
         constraints(max_number_constraints),
         objective_gradient(max_number_variables),
         constraint_jacobian(max_number_constraints, max_number_variables) {
   }
};

#endif // UNO_EVALUATIONS_H
