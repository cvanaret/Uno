// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ElasticVariables.hpp"
#include "model/Model.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   ElasticVariables::ElasticVariables(size_t number_positive_variables, size_t number_negative_variables):
         positive(number_positive_variables), negative(number_negative_variables) { }

   ElasticVariables ElasticVariables::generate(const Model& model) {
      const ElasticVariablesSizes sizes = ElasticVariables::count(model);
      ElasticVariables elastic_variables(sizes.positive, sizes.negative);

      // generate elastic variables to relax the constraints
      size_t elastic_index = model.number_variables;
      for (size_t constraint_index: Range(model.number_constraints)) {
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

   ElasticVariablesSizes ElasticVariables::count(const Model& model) {
      ElasticVariablesSizes number_elastic_variables;
      // if the subproblem uses slack variables, the bounds of the constraints are [0, 0]
      for (size_t constraint_index: Range(model.number_constraints)) {
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