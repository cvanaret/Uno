#include "SubproblemFactory.hpp"
#include "QPSubproblem.hpp"
#include "LPSubproblem.hpp"
#include "BarrierSubproblem.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(const Problem& problem, size_t max_number_variables, const Options& options) {
   const std::vector<std::string> possible_methods = {"QP", "LP", "barrier"};
   const std::string subproblem_type = options.at("subproblem");
   // active-set methods
   if (subproblem_type == "QP") {
      return std::make_unique<QPSubproblem>(problem, max_number_variables, options);
   }
   else if (subproblem_type == "LP") {
      return std::make_unique<LPSubproblem>(problem, max_number_variables, options);
   }
   // interior point method
   else if (subproblem_type == "barrier") {
      return std::make_unique<BarrierSubproblem>(problem, max_number_variables, options);
   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " is not supported");
}