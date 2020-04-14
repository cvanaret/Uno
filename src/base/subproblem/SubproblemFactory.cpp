#include "SubproblemFactory.hpp"
#include "SQP.hpp"
#include "SLP.hpp"
#include "SLPEQP.hpp"
#include "InteriorPoint.hpp"
#include "QPSolverFactory.hpp"

std::shared_ptr<Subproblem> SubproblemFactory::create(const std::string& type, QPSolver& solver, std::map<std::string, std::string> default_values) {
    if (type == "SQP") {
        return std::make_shared<SQP>(solver);
    }
    else if (type == "SLP") {
        return std::make_shared<SLP>(solver);
    }
    else if (type == "SLPEQP") {
        return std::make_shared<SLPEQP>(solver);
    }
    else if (type == "IPM") {
        return std::make_shared<InteriorPoint>();
    }
    //else if (type == "AL") {
    //    return std::make_shared<AugmentedLagrangian>();
    //}
    else {
        throw std::invalid_argument("LocalApproximation type " + type + " does not exist");
    }
}
