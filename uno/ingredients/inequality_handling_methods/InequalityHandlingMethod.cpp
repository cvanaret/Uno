// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const std::string& hessian_model, const std::string& regularization_strategy,
      size_t number_variables, const Options& options) :
         hessian_model(HessianModelFactory::create(hessian_model)),
         regularization_strategy(RegularizationStrategyFactory::create(regularization_strategy, options)),
         initial_point(number_variables) {
   }

   size_t InequalityHandlingMethod::get_hessian_evaluation_count() const {
      return this->hessian_model->evaluation_count;
   }
} // namespace