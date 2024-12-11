// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InequalityHandlingMethod.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const std::string& hessian_model, size_t dimension, size_t number_hessian_nonzeros, bool convexify,
         const Options& options) :
         hessian_model(HessianModelFactory::create(hessian_model, dimension, number_hessian_nonzeros, convexify, options)) {
   }

   void InequalityHandlingMethod::set_trust_region_radius(double new_trust_region_radius) {
      assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
      this->trust_region_radius = new_trust_region_radius;
   }

   size_t InequalityHandlingMethod::get_hessian_evaluation_count() const {
      return this->hessian_model->evaluation_count;
   }
} // namespace