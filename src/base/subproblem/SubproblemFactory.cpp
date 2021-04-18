#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "Sl1QP.hpp"
#include "SLPEQP.hpp"
#include "InteriorPoint.hpp"
#include "QPSolverFactory.hpp"

std::unique_ptr<Subproblem> SubproblemFactory::create(Problem& problem, const std::string& type, std::map<std::string, std::string> options, bool use_trust_region, bool scale_residuals) {
    std::vector<std::string> possible_methods = {"SQP", "SLP", "Sl1QP", "SLPEQP", "IPM"};
    /* active-set methods */
    if (type == "SQP") {
        return std::make_unique<SQP>(problem, options["QP_solver"], options["hessian"], use_trust_region, scale_residuals);
    }
    else if (type == "SLP") {
        return std::make_unique<SLP>(problem, options["QP_solver"], use_trust_region, scale_residuals);
    }
    else if (type == "Sl1QP") {
        double initial_parameter = std::stod(options["Sl1QP_initial_parameter"]);
        return std::make_unique<Sl1QP>(problem, options["QP_solver"], options["hessian"], use_trust_region, scale_residuals, initial_parameter);
    }
    else if (type == "SLPEQP") {
        return std::make_unique<SLPEQP>(problem, options["QP_solver"], options["linear_solver"], options["hessian"], use_trust_region, scale_residuals);
    }
//    else if (type == "SLPEQP") {
//          if (use_trust_region) {
//             return std::make_unique<SLPEQP_TR>(problem, options["LP_solver"], options["hessian"], use_trust_region, scale_residuals);
//          }
//          else {
//             return std::make_unique<SLPEQP_l2>(problem, options["hessian"], use_trust_region, scale_residuals);
//          }
//    }
    /* interior point method */
    else if (type == "IPM") {
        return std::make_unique<InteriorPoint>(problem, options["linear_solver"], options["hessian"], use_trust_region, scale_residuals);
    }
    else {
        throw std::invalid_argument("Subproblem method " + type + " does not exist. The possible options are: " + join(possible_methods, ", "));
    }
}
