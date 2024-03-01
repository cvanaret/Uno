// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ModelFactory.hpp"
#include "HomogeneousEqualityConstrainedModel.hpp"
#include "ScaledModel.hpp"
#include "BoundRelaxedModel.hpp"
#include "preprocessing/Scaling.hpp"

// note: transfer ownership of the pointer
std::unique_ptr<Model> ModelFactory::reformulate(std::unique_ptr<Model> model, Iterate& initial_iterate, const Options& options) {
   // optional: scale the problem using the evaluations at the first iterate
   if (options.get_string("scale_functions") == "yes") {
      model = std::make_unique<ScaledModel>(std::move(model), initial_iterate, options);
   }

   // if an equality-constrained problem is required (e.g. interior points or AL), reformulate the model with slacks
   if (options.get_string("subproblem") == "primal_dual_interior_point") {
      // generate an equality-constrained model by:
      // - introducing slacks in inequality constraints
      // - subtracting the (possibly nonzero) RHS of equality constraints
      model = std::make_unique<HomogeneousEqualityConstrainedModel>(std::move(model));
      // slightly relax the bound constraints
      model = std::make_unique<BoundRelaxedModel>(std::move(model), options);
      initial_iterate.set_number_variables(model->number_variables);
   }
   return model;
}
