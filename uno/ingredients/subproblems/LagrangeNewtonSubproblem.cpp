// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "LagrangeNewtonSubproblem.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "tools/Options.hpp"

namespace uno {
   LagrangeNewtonSubproblem::LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, bool use_regularization, double trust_region_radius, const Options& options):
      problem(problem),
      current_iterate(current_iterate),
      current_multipliers(current_multipliers),
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), problem.number_variables,
            problem.number_hessian_nonzeros() + (use_regularization ? problem.number_variables : 0), use_regularization, options)),
      trust_region_radius(trust_region_radius) { }

   LagrangeNewtonSubproblem::~LagrangeNewtonSubproblem() { }

   const SymmetricMatrix<size_t, double>& LagrangeNewtonSubproblem::get_lagrangian_hessian() const {
      return this->hessian_model->hessian;
   }
} // namespace