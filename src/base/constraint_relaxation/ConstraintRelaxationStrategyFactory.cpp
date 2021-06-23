#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const std::string& constraint_relaxation_type,
      Problem& problem, Subproblem& subproblem, const std::map<std::string, std::string>& options) {
   if (constraint_relaxation_type == "feasibility-restoration") {
      return std::make_unique<FeasibilityRestoration>(subproblem, options);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      return std::make_unique<l1Relaxation>(problem, subproblem, options);
   }
   else {
      throw std::invalid_argument("ConstraintRelaxationStrategy type " + constraint_relaxation_type + " does not exist");
   }
}