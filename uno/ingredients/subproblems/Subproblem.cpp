// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Subproblem.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"

namespace uno {
   Subproblem::Subproblem(const std::string& hessian_model, size_t dimension, size_t number_hessian_nonzeros, bool convexify,
         const Options& options) :
         hessian_model(HessianModelFactory::create(hessian_model, dimension, number_hessian_nonzeros, convexify, options)) {
   }

   void Subproblem::set_trust_region_radius(double new_trust_region_radius) {
      assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
      this->trust_region_radius = new_trust_region_radius;
   }

   const SymmetricMatrix<size_t, double>& Subproblem::get_lagrangian_hessian() const {
      return this->hessian_model->hessian;
   }

   size_t Subproblem::get_hessian_evaluation_count() const {
      return this->hessian_model->evaluation_count;
   }
} // namespace