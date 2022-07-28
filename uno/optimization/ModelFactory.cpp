// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "EqualityConstrainedModel.hpp"

std::unique_ptr<Model> ModelFactory::create(const std::string& problem_name, const Options& options) {
   // create the AMPL model
   std::unique_ptr<Model> model = std::make_unique<AMPLModel>(problem_name);

   // if an equality-constrained problem is required (e.g. barrier or AL), reformulate the model
   if (options.get_string("subproblem") == "barrier") {
      // transfer ownership of the pointer
      model = std::make_unique<EqualityConstrainedModel>(std::move(model));
   }
   return model;
}