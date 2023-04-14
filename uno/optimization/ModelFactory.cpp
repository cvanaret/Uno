// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "EqualityConstrainedModel.hpp"
#include "ScaledModel.hpp"
#include "BoundRelaxedModel.hpp"
#include "preprocessing/Scaling.hpp"

// note: transfer ownership of the pointer
std::unique_ptr<Model> ModelFactory::reformulate(std::unique_ptr<Model> model, Iterate& first_iterate, const Options& options) {
   // optional: scale the problem using the evaluations at the first iterate
   if (options.get_string("scale_functions") == "yes") {
      model = std::make_unique<ScaledModel>(std::move(model), first_iterate, options);
   }

   // if an equality-constrained problem is required (e.g. barrier or AL), reformulate the model with slacks
   if (options.get_string("subproblem") == "barrier") {
      if (not model->inequality_constraints.empty()) {
         // introduce slacks to obtain equality constraints
         model = std::make_unique<EqualityConstrainedModel>(std::move(model));
      }
      // slightly relax the bound constraints
      model = std::make_unique<BoundRelaxedModel>(std::move(model), options);
      first_iterate.set_number_variables(model->number_variables);
   }
   return model;
}