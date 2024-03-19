// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "HomogeneousEqualityConstrainedModel.hpp"
#include "ScaledModel.hpp"
#include "BoundRelaxedModel.hpp"
#include "preprocessing/Scaling.hpp"

// note: ownership of the pointer is transferred
std::unique_ptr<Model> ModelFactory::reformulate(std::unique_ptr<Model> model, Iterate& initial_iterate, const Options& options) {
   // optional: scale the problem using the evaluations at the first iterate
   if (options.get_string("scale_functions") == "yes") {
      model = std::make_unique<ScaledModel>(std::move(model), initial_iterate, options);
   }

   if (options.get_string("subproblem") == "primal_dual_interior_point") {
      // if an equality-constrained problem is required (e.g. interior points or AL), reformulate the model with slacks
      model = std::make_unique<HomogeneousEqualityConstrainedModel>(std::move(model));
      // slightly relax the bound constraints
      model = std::make_unique<BoundRelaxedModel>(std::move(model), options);
      // add the slacks to the initial iterate
      initial_iterate.set_number_variables(model->number_variables);
   }
   return model;
}
