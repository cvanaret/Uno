// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONS_H
#define UNO_EVALUATIONS_H

#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // forward declaration
   class Model;

   struct Evaluations {
      double objective{INF<double>}; /*!< Objective value */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      Vector<double> objective_gradient; /*!< Sparse Jacobian of the objective */

      Evaluations(size_t number_variables, size_t number_constraints);

      void evaluate_objective(const Model& model, const Vector<double>& primals);
      void evaluate_constraints(const Model& model, const Vector<double>& primals);
      void evaluate_objective_gradient(const Model& model, const Vector<double>& primals);

      void compute_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const;
   };
} // namespace

#endif // UNO_EVALUATIONS_H