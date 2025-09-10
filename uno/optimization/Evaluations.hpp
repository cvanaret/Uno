// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONS_H
#define UNO_EVALUATIONS_H

#include <vector>
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
      struct Evaluations {
            double objective{INF<double>}; /*!< Objective value */
            Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
            std::vector<double> linearized_constraints;
            Vector<double> objective_gradient; /*!< Sparse Jacobian of the objective */

            Evaluations(size_t number_variables, size_t number_constraints):
                  constraints(number_constraints),
                  linearized_constraints(number_constraints),
                  objective_gradient(number_variables) {
            }
      };
} // namespace

#endif // UNO_EVALUATIONS_H