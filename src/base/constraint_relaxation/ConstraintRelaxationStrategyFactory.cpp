#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const std::string& constraint_relaxation_type,
      Problem& problem, const std::map<std::string, std::string>& options, bool use_trust_region) {
   if (constraint_relaxation_type == "feasibility-restoration") {
      return std::make_unique<FeasibilityRestoration>(problem, options, use_trust_region);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      return std::make_unique<l1Relaxation>(problem, options, use_trust_region);
   }
   else {
      throw std::invalid_argument("ConstraintRelaxationStrategy type " + constraint_relaxation_type + " does not exist");
   }
}