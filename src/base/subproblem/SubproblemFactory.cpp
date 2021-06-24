#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "InteriorPoint.hpp"
#include "QPSolverFactory.hpp"
#include "Vector.hpp"

std::unique_ptr<Subproblem>
SubproblemFactory::create(const Problem& problem, size_t number_variables, const std::string& subproblem_type, const std::map<std::string,
      std::string>& options, bool use_trust_region) {
   const std::vector<std::string> possible_methods = {"SQP", "SLP", "IPM"};
   /* active-set methods */
   if (subproblem_type == "SQP") {
      return std::make_unique<SQP>(problem, number_variables, problem.number_constraints, options.at("QP_solver"), options.at("hessian"),
            use_trust_region);
   }
   else if (subproblem_type == "SLP") {
      return std::make_unique<SLP>(number_variables, problem.number_constraints, options.at("QP_solver"));
   }
//    else if (type == "SLPEQP") {
//          if (use_trust_region) {
//             return std::make_unique<SLPEQP_TR>(problem, options.at("LP_solver"], options.at("hessian"], use_trust_region, scale_residuals);
//          }
//          else {
//             return std::make_unique<SLPEQP_l2>(problem, options.at("hessian"], use_trust_region, scale_residuals);
//          }
//    }
   /* interior point method */
   else if (subproblem_type == "IPM") {
      return std::make_unique<InteriorPoint>(problem, number_variables, problem.number_constraints, options.at("linear_solver"), options.at
      ("hessian"), use_trust_region);
   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " does not exist.");
}
