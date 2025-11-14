// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Scaling.hpp"
#include "Iterate.hpp"
#include "linear_algebra/Norm.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "tools/Logger.hpp"

namespace uno {
   Scaling::Scaling(const Iterate& initial_iterate, const EvaluationSpace& evaluation_space, double gradient_threshold):
         gradient_threshold(gradient_threshold), objective_scaling(1.),
         constraint_scaling(initial_iterate.evaluations.constraints.size(), 1.) {
      // objective
      this->objective_scaling = std::min(1., this->gradient_threshold / norm_inf(initial_iterate.evaluations.objective_gradient));

      // constraints
      Vector<double> row_norms(initial_iterate.evaluations.constraints.size());
      evaluation_space.compute_constraint_jacobian_norms(row_norms);
      for (size_t constraint_index: Range(this->constraint_scaling.size())) {
         this->constraint_scaling[constraint_index] = std::min(1., this->gradient_threshold / row_norms[constraint_index]);
      }

      DEBUG << "Objective scaling: " << this->objective_scaling << '\n';
      DEBUG << "Constraint scaling: "; print_vector(DEBUG, this->constraint_scaling);
   }

   double Scaling::get_objective_scaling() const {
      return this->objective_scaling;
   }

   double Scaling::get_constraint_scaling(size_t constraint_index) const {
      assert(constraint_index < this->constraint_scaling.size() && "The constraint index is not valid.");
      return this->constraint_scaling[constraint_index];
   }
} // namespace