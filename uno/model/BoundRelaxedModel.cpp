// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include "BoundRelaxedModel.hpp"
#include "symbolic/Collection.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   BoundRelaxedModel::BoundRelaxedModel(const Model& original_model, const Options& options):
         Model(original_model.name + " -> bounds relaxed", original_model.number_variables, original_model.number_constraints,
            original_model.optimization_sense, original_model.lagrangian_sign_convention, original_model.base_indexing),
         model(original_model),
         relaxation_factor(options.get_double("bound_relaxation_factor")),
         constraint_violation_tolerance(options.get_double("constraint_violation_tolerance")),
         relaxed_variables_lower_bounds(this->model.get_variables_lower_bounds()),
         relaxed_variables_upper_bounds(this->model.get_variables_upper_bounds()),
         relaxed_constraints_lower_bounds(this->model.get_constraints_lower_bounds()),
         relaxed_constraints_upper_bounds(this->model.get_constraints_upper_bounds()) {
      if (this->relaxation_factor <= 0.) {
         return; // no relaxation unless bound_relax_factor > 0
      }

      const auto relax_lower = [&](double& bound) {
         if (is_finite(bound)) {
            bound -= bound_relaxation(bound);
         }
      };
      const auto relax_upper = [&](double& bound) {
         if (is_finite(bound)) {
            bound += bound_relaxation(bound);
         }
      };

      // relaxed variables bounds
      for (size_t variable_index: Range(this->number_variables)) {
         relax_lower(this->relaxed_variables_lower_bounds[variable_index]);
         relax_upper(this->relaxed_variables_upper_bounds[variable_index]);
      }

      // relaxed constraints bounds
      for (size_t constraint_index: this->model.get_inequality_constraints()) {
         relax_lower(this->relaxed_constraints_lower_bounds[constraint_index]);
         relax_upper(this->relaxed_constraints_upper_bounds[constraint_index]);
      }
   }

   // relaxation amount, as in Ipopt's OrigIpoptNLP::relax_bounds:
   // min(constr_viol_tol, |bound_relax_factor| * max(1, |bound|))
   double BoundRelaxedModel::bound_relaxation(double bound) const {
      return std::min(this->constraint_violation_tolerance, std::abs(this->relaxation_factor) * std::max(1., std::abs(bound)));
   }
} // namespace