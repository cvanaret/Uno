#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "InteriorPoint.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(const Problem& problem, size_t max_number_variables, const Options& options) {
   const std::vector<std::string> possible_methods = {"SQP", "SLP", "IPM"};
   const std::string subproblem_type = options.at("subproblem");
   // active-set methods
   if (subproblem_type == "SQP") {
      return std::make_unique<SQP>(problem, max_number_variables, options);
   }
   else if (subproblem_type == "SLP") {
      return std::make_unique<SLP>(problem, max_number_variables, options);
   }
   // interior point method
   else if (subproblem_type == "IPM") {
      return std::make_unique<InteriorPoint>(problem, max_number_variables, options);
   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " is not supported");
}
