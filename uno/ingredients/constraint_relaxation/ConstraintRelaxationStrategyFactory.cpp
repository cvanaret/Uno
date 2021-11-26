#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const std::string& constraint_relaxation_type,
      Problem& problem, const Options& options) {
   if (constraint_relaxation_type == "feasibility-restoration") {
      return std::make_unique<FeasibilityRestoration>(problem, options);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      const double initial_parameter = stod(options.at("l1_relaxation_initial_parameter"));
      const double decrease_factor = stod(options.at("l1_relaxation_decrease_factor"));
      const double epsilon1 = stod(options.at("l1_relaxation_epsilon1"));
      const double epsilon2 = stod(options.at("l1_relaxation_epsilon2"));
      l1RelaxationParameters parameters({initial_parameter, decrease_factor, epsilon1, epsilon2});
      return std::make_unique<l1Relaxation>(problem, parameters, options);
   }
   throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
}