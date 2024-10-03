// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Scaling.hpp"
#include "linear_algebra/Norm.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   Scaling::Scaling(size_t number_constraints, double gradient_threshold):
         gradient_threshold(gradient_threshold), objective_scaling(1.), constraint_scaling(number_constraints, 1.) {
   }

   void Scaling::compute(const SparseVector<double>& objective_gradient, const RectangularMatrix<double>& constraint_jacobian) {
      // objective
      this->objective_scaling = std::min(1., this->gradient_threshold / norm_inf(objective_gradient));

      // constraints
      for (size_t constraint_index: Range(this->constraint_scaling.size())) {
         this->constraint_scaling[constraint_index] = std::min(1., this->gradient_threshold / norm_inf(constraint_jacobian[constraint_index]));
      }
      DEBUG2 << "Objective scaling: " << this->objective_scaling << '\n';
      DEBUG2 << "Constraint scaling: "; print_vector(DEBUG2, this->constraint_scaling);
   }

   double Scaling::get_objective_scaling() const {
      return this->objective_scaling;
   }

   double Scaling::get_constraint_scaling(size_t constraint_index) const {
      assert(constraint_index < this->constraint_scaling.size() && "The constraint index is not valid.");
      return this->constraint_scaling[constraint_index];
   }
} // namespace
