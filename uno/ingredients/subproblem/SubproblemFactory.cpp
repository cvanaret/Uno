#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "InteriorPoint.hpp"

std::unique_ptr <Subproblem> SubproblemFactory::create(const Problem& problem, size_t number_variables, const std::string& subproblem_type,
      const Options& options, bool use_trust_region) {
   const std::vector <std::string> possible_methods = {"SQP", "SLP", "IPM"};
   /* active-set methods */
   if (subproblem_type == "SQP") {
      const std::string& QP_solver_name = options.at("QP_solver");
      // determine the sparse matrix format
      if (QP_solver_name == "BQPD") {
         return std::make_unique<SQP<BQPDSolver>>(problem, number_variables, problem.number_constraints, options.at("hessian"), use_trust_region);
      }
      else {
         assert(false && "SubproblemFactory::create: unknown QP solver");
      }
   }
   else if (subproblem_type == "SLP") {
      const std::string& LP_solver_name = options.at("LP_solver");
      if (LP_solver_name == "BQPD") {
         return std::make_unique<SLP<BQPDSolver>>(problem, number_variables, problem.number_constraints);
      }
      else {
         assert(false && "SubproblemFactory::create: unknown LP solver");
      }
   }
      /* interior point method */
   else if (subproblem_type == "IPM") {
      const std::string& linear_solver_name = options.at("linear_solver");
      const double initial_barrier_parameter = std::stod(options.at("initial_barrier_parameter"));
      const double default_multiplier = std::stod(options.at("default_multiplier"));
      const double tolerance = std::stod(options.at("tolerance"));
      // determine the sparse matrix format
      if (linear_solver_name == "MA57") {
         return std::make_unique<InteriorPoint<MA57Solver>>(problem, number_variables, problem.number_constraints, options.at("hessian"),
               initial_barrier_parameter, default_multiplier, tolerance, use_trust_region);
      }
      else if (linear_solver_name == "PARDISO") {
         return std::make_unique<InteriorPoint<PardisoSolver>>(problem, number_variables, problem.number_constraints, options.at("hessian"),
               initial_barrier_parameter, default_multiplier, tolerance, use_trust_region);
      }
      else {
         assert(false && "SubproblemFactory::create: unknown linear solver");
      }

   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " does not exist.");
}
