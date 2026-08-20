// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Scaling.hpp"
#include "Model.hpp"
#include "linear_algebra/Norm.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

namespace uno {
   Scaling::Scaling(size_t number_constraints, double gradient_threshold):
         gradient_threshold(gradient_threshold), objective_scaling(1.), constraint_scaling(number_constraints, 1.) {
   }

   void Scaling::compute(const Model& model, const Vector<double>& objective_gradient, const Vector<double>& jacobian_values) {
      // objective
      this->objective_scaling = std::min(1., this->gradient_threshold / norm_inf(objective_gradient));

      // constraints
      // compute the inf norm of each row of the Jacobian
      Vector<double> norm_inf_constraints(model.number_constraints, 0.);
      const auto& jacobian_row_indices = model.get_jacobian_row_indices();
      for (size_t nonzero_index: Range(model.number_jacobian_nonzeros())) {
         const size_t constraint_index = static_cast<size_t>(jacobian_row_indices[nonzero_index]);
         norm_inf_constraints[constraint_index] = std::max(norm_inf_constraints[constraint_index], std::abs(jacobian_values[nonzero_index]));
      }
      for (size_t constraint_index: Range(model.number_constraints)) {
         const double row_norm = norm_inf_constraints[constraint_index];
         this->constraint_scaling[constraint_index] = (row_norm > 0.) ? std::min(1., this->gradient_threshold / row_norm) : 1.;
      }
      DEBUG2 << "Objective scaling: " << this->objective_scaling << '\n';
      DEBUG2 << "Constraint scaling: " << view(this->constraint_scaling) << '\n';
   }

   double Scaling::get_objective_scaling() const {
      return this->objective_scaling;
   }

   double Scaling::get_constraint_scaling(size_t constraint_index) const {
      assert(constraint_index < this->constraint_scaling.size() && "The constraint index is not valid.");
      return this->constraint_scaling[constraint_index];
   }
} // namespace