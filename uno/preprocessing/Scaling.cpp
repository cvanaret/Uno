// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Scaling.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"

Scaling::Scaling(size_t number_constraints, double gradient_threshold):
      gradient_threshold(gradient_threshold),
      objective_scaling(1.),
      constraint_scaling(number_constraints, 1.) {
}

void Scaling::compute(const SparseVector<double>& objective_gradient, const RectangularMatrix<double>& constraint_jacobian) {
   // set the objective scaling
   this->objective_scaling = std::min(1., this->gradient_threshold / norm_inf(objective_gradient));

   // set the constraints scaling
   for (size_t j: Range(this->constraint_scaling.size())) {
      this->constraint_scaling[j] = std::min(1., this->gradient_threshold / norm_inf(constraint_jacobian[j]));
   }
   DEBUG2 << "Objective scaling: " << this->objective_scaling << '\n';
   DEBUG2 << "Constraint scaling: "; print_vector(DEBUG2, this->constraint_scaling);
}

double Scaling::get_objective_scaling() const {
   return this->objective_scaling;
}

double Scaling::get_constraint_scaling(size_t j) const {
   return this->constraint_scaling[j];
}