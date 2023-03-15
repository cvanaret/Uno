// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "EqualityConstrainedModel.hpp"
#include "ScaledModel.hpp"
#include "BoundRelaxedModel.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "preprocessing/Scaling.hpp"

std::unique_ptr<Model> ModelFactory::reformulate(const AMPLModel& ampl_model, Iterate& first_iterate, const Options& options) {
   std::unique_ptr<Model> model = std::make_unique<ScaledModel>(ampl_model, first_iterate, options);

   // if an equality-constrained problem is required (e.g. barrier or AL), reformulate the model
   if (options.get_string("subproblem") == "barrier") {
      // transfer ownership of the pointers
      // introduce slacks to obtain equality constraints
      model = std::make_unique<EqualityConstrainedModel>(std::move(model));
      // slightly relax the bound constraints
      model = std::make_unique<BoundRelaxedModel>(std::move(model), options);
      first_iterate.set_number_variables(model->number_variables);
   }
   return model;
}