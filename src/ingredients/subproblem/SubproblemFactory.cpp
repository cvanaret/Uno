#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "InteriorPoint.hpp"
#include "QPSolverFactory.hpp"
#include "Vector.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(const Problem& problem, size_t number_variables, const std::string& subproblem_type,
      const Options& options, bool use_trust_region) {
   const std::vector<std::string> possible_methods = {"SQP", "SLP", "IPM"};
   /* active-set methods */
   if (subproblem_type == "SQP") {
      const std::string& QP_solver_name = options.at("QP_solver");
      // determine the sparse matrix format
      if (QP_solver_name == "BQPD") {
         //std::cout << "Use CSC\n";
      }
      else {
         assert(false && "SubproblemFactory::create: unknown QP solver");
      }
      return std::make_unique<SQP>(problem, number_variables, problem.number_constraints, QP_solver_name, options.at("hessian"),
            use_trust_region);
   }
   else if (subproblem_type == "SLP") {
      const std::string& QP_solver_name = options.at("QP_solver");
      return std::make_unique<SLP>(number_variables, problem.number_constraints, QP_solver_name);
   }
   /* interior point method */
   else if (subproblem_type == "IPM") {
      const std::string& linear_solver_name = options.at("linear_solver");
      // determine the sparse matrix format
      if (linear_solver_name == "MA57") {
         return std::make_unique<InteriorPoint<COOSymmetricMatrix> >(problem, number_variables, problem.number_constraints, linear_solver_name,
               options.at("hessian"), use_trust_region);
      }
      else if (linear_solver_name == "PARDISO") {
         return std::make_unique<InteriorPoint<CSCSymmetricMatrix> >(problem, number_variables, problem.number_constraints, linear_solver_name,
               options.at("hessian"), use_trust_region);
      }
      else {
         assert(false && "SubproblemFactory::create: unknown QP solver");
      }

   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " does not exist.");
}
