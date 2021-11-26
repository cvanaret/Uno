#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "InteriorPoint.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(const Problem& problem, size_t max_number_variables, const Options& options) {
   const std::vector<std::string> possible_methods = {"SQP", "SLP", "IPM"};
   const std::string subproblem_type = options.at("subproblem");
   const std::string mechanism_type = options.at("mechanism");
   const bool use_trust_region = (mechanism_type == "TR");

   // active-set methods
   if (subproblem_type == "SQP") {
      const std::string& hessian_model = options.at("hessian_model");
      const std::string& QP_solver_name = options.at("QP_solver");
      const std::string& sparse_format = options.at("sparse_format");
      return std::make_unique<SQP>(problem, max_number_variables, problem.number_constraints, hessian_model, QP_solver_name, sparse_format,
            use_trust_region);
   }
   else if (subproblem_type == "SLP") {
      const std::string& LP_solver_name = options.at("LP_solver");
      return std::make_unique<SLP>(problem, max_number_variables, problem.number_constraints, LP_solver_name);
   }
   // interior point method
   else if (subproblem_type == "IPM") {
      const std::string& hessian_model = options.at("hessian_model");
      const std::string& linear_solver_name = options.at("linear_solver");
      const std::string& sparse_format = options.at("sparse_format");
      const double initial_barrier_parameter = std::stod(options.at("initial_barrier_parameter"));
      const double default_multiplier = std::stod(options.at("default_multiplier"));
      const double tolerance = std::stod(options.at("tolerance"));
      return std::make_unique<InteriorPoint>(problem, max_number_variables, problem.number_constraints, hessian_model, linear_solver_name,
            sparse_format, initial_barrier_parameter, default_multiplier, tolerance, use_trust_region);
   }
   throw std::invalid_argument("Subproblem method " + subproblem_type + " is not supported");
}
