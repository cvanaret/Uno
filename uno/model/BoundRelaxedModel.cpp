// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "BoundRelaxedModel.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace uno {
   BoundRelaxedModel::BoundRelaxedModel(const Model& original_model, const Options& options):
         Model(original_model.name + " -> bounds relaxed", original_model.number_variables, original_model.number_constraints,
            original_model.optimization_sense),
         model(original_model),
         relaxation_factor(options.get_double("primal_tolerance")) {
   }

   double BoundRelaxedModel::variable_lower_bound(size_t variable_index) const {
      const double lower_bound = this->model.variable_lower_bound(variable_index);
      return lower_bound - this->relaxation_factor * std::max(1., std::abs(lower_bound));
   }

   double BoundRelaxedModel::variable_upper_bound(size_t variable_index) const {
      const double upper_bound = this->model.variable_upper_bound(variable_index);
      return upper_bound + this->relaxation_factor * std::max(1., std::abs(upper_bound));
   }
} // namespace