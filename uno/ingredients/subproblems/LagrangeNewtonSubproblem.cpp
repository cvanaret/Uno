// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "LagrangeNewtonSubproblem.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "tools/Options.hpp"

namespace uno {
   LagrangeNewtonSubproblem::LagrangeNewtonSubproblem(const OptimizationProblem& problem, const Iterate& current_iterate, size_t number_variables,
         size_t number_hessian_nonzeros, bool use_regularization, double trust_region_radius, const Options& options):
         problem(problem),
         current_iterate(current_iterate),
         hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), number_variables,
               number_hessian_nonzeros + (use_regularization ? number_variables : 0), use_regularization, options)),
         trust_region_radius(trust_region_radius) { }

   LagrangeNewtonSubproblem::~LagrangeNewtonSubproblem() { }

   void LagrangeNewtonSubproblem::evaluate_functions(SymmetricMatrix<size_t, double>& hessian, const WarmstartInformation& warmstart_information) {
      // TODO evaluate functions
   }
} // namespace