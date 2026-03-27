// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "BoundRelaxedModel.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace uno {
   BoundRelaxedModel::BoundRelaxedModel(const Model& original_model, const Options& options):
         Model(original_model.name + " -> bounds relaxed", original_model.number_variables, original_model.number_constraints,
            original_model.optimization_sense, original_model.lagrangian_sign_convention),
         model(original_model),
         relaxation_factor(options.get_double("primal_tolerance")),
         relaxed_variables_lower_bounds(this->number_variables),
         relaxed_variables_upper_bounds(this->number_variables) {
      // compute the relaxed bounds
      const auto& variables_lower_bounds = this->model.get_variables_lower_bounds();
      const auto& variables_upper_bounds = this->model.get_variables_upper_bounds();
      for (size_t variable_index: Range(this->number_variables)) {
         const double lower_bound = variables_lower_bounds[variable_index];
         this->relaxed_variables_lower_bounds[variable_index] = lower_bound - this->relaxation_factor *
            std::max(1., std::abs(lower_bound));
         const double upper_bound = variables_upper_bounds[variable_index];
         this->relaxed_variables_upper_bounds[variable_index] = upper_bound + this->relaxation_factor *
            std::max(1., std::abs(upper_bound));;
      }
   }


} // namespace