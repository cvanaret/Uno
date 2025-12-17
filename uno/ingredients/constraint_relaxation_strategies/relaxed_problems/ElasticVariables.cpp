// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ElasticVariables.hpp"
#include "model/Model.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   ElasticVariables::ElasticVariables(size_t number_positive_variables, size_t number_negative_variables):
         positive(number_positive_variables), negative(number_negative_variables) { }

   ElasticVariables ElasticVariables::generate(const Model& model, bool relax_linear_constraints) {
      const ElasticVariablesSizes sizes = ElasticVariables::count(model, relax_linear_constraints);
      ElasticVariables elastic_variables(sizes.positive, sizes.negative);

      // generate elastic variables to relax the constraints
      size_t elastic_index = model.number_variables;
      const auto& constraints = relax_linear_constraints ?
         static_cast<const Collection<size_t>&>(Range(model.number_constraints)) : model.get_nonlinear_constraints();
      for (size_t constraint_index: constraints) {
         if (is_finite(model.constraint_upper_bound(constraint_index))) {
            // nonnegative variable p that captures the positive part of the constraint violation
            elastic_variables.positive.insert(constraint_index, elastic_index);
            elastic_index++;
         }
         if (is_finite(model.constraint_lower_bound(constraint_index))) {
            // nonpositive variable n that captures the negative part of the constraint violation
            elastic_variables.negative.insert(constraint_index, elastic_index);
            elastic_index++;
         }
      }
      return elastic_variables;
   }

   ElasticVariablesSizes ElasticVariables::count(const Model& model, bool relax_linear_constraints) {
      ElasticVariablesSizes number_elastic_variables{};
      const auto& constraints = relax_linear_constraints ?
         static_cast<const Collection<size_t>&>(Range(model.number_constraints)) : model.get_nonlinear_constraints();
      for (size_t constraint_index: constraints) {
         if (is_finite(model.constraint_upper_bound(constraint_index))) {
            number_elastic_variables.positive++;
         }
         if (is_finite(model.constraint_lower_bound(constraint_index))) {
            number_elastic_variables.negative++;
         }
      }
      return number_elastic_variables;
   }
} // namespace