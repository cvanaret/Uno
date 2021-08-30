#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

size_t ConstraintRelaxationStrategyFactory::get_number_variables(std::string_view constraint_relaxation_type, const Problem& problem) {
   if (constraint_relaxation_type == "feasibility-restoration") {
      return FeasibilityRestoration::get_number_variables(problem);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      return l1Relaxation::get_number_variables(problem);
   }
   else {
      throw std::invalid_argument("ConstraintRelaxationStrategy type does not exist");
   }
}

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(std::string_view constraint_relaxation_type,
      Problem& problem, Subproblem& subproblem, const Options& options) {
   if (constraint_relaxation_type == "feasibility-restoration") {
      return std::make_unique<FeasibilityRestoration>(subproblem, options);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      return std::make_unique<l1Relaxation>(problem, subproblem, options);
   }
   else {
      throw std::invalid_argument("ConstraintRelaxationStrategy type does not exist");
   }
}