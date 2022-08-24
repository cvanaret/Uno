// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SubproblemFactory.hpp"
#include "ingredients/subproblem/active_set/QPSubproblem.hpp"
#include "ingredients/subproblem/active_set/LPSubproblem.hpp"
#include "ingredients/subproblem/interior_point/InfeasibleInteriorPointSubproblem.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros,
      const Options& options) {
   const std::vector<std::string> possible_methods = {"QP", "LP", "barrier"};
   const std::string subproblem_type = options.get_string("subproblem");
   // active-set methods
   if (subproblem_type == "QP") {
      return std::make_unique<QPSubproblem>(max_number_variables, max_number_constraints, max_number_hessian_nonzeros, options);
   }
   else if (subproblem_type == "LP") {
      return std::make_unique<LPSubproblem>(max_number_variables, max_number_constraints, options);
   }
   // interior-point method
   else if (subproblem_type == "barrier") {
      return std::make_unique<InfeasibleInteriorPointSubproblem>(max_number_variables, max_number_constraints, max_number_hessian_nonzeros, options);
   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " is not supported");
}