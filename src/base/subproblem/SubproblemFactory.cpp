#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "Sl1QP.hpp"
//#include "SLPEQP.hpp"
#include "InteriorPoint.hpp"
#include "QPSolverFactory.hpp"

std::shared_ptr<Subproblem> SubproblemFactory::create(Problem& problem, const std::string& type, std::map<std::string, std::string> default_values, bool use_trust_region, bool scale_residuals) {
    /* active-set methods */
    if (type == "SQP") {
        return std::make_shared<SQP>(problem, default_values["QP_solver"], default_values["hessian"], use_trust_region, scale_residuals);
    }
    else if (type == "SLP") {
        return std::make_shared<SLP>(problem, default_values["QP_solver"], use_trust_region, scale_residuals);
    }
    else if (type == "Sl1QP") {
        return std::make_shared<Sl1QP>(problem, default_values["QP_solver"], default_values["hessian"], use_trust_region, scale_residuals);
    }
//    else if (type == "SLPEQP") {
//        return std::make_shared<SLPEQP>(problem, default_values["QP_solver"], default_values["hessian"]);
//    }
    /* interior point method */
    else if (type == "IPM") {
        return std::make_shared<InteriorPoint>(problem, default_values["hessian"], use_trust_region, scale_residuals);
    }
    else {
        throw std::invalid_argument("Subproblem method " + type + " does not exist");
    }
}
