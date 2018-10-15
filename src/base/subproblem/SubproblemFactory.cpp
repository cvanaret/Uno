#include "SubproblemFactory.hpp"
#include "QPApproximation.hpp"
#include "InteriorPoint.hpp"
#include "AugmentedLagrangian.hpp"
#include "QPSolverFactory.hpp"

std::shared_ptr<Subproblem> SubproblemFactory::create(const std::string& type, QPSolver& solver, std::map<std::string, std::string> default_values) {
    if (type == "QP") {
        return std::make_shared<QPApproximation>(solver);
    }
    else if (type == "IPM") {
        return std::make_shared<InteriorPoint>();
    }
    else if (type == "AL") {
        return std::make_shared<AugmentedLagrangian>();
    }
    else {
        throw std::invalid_argument("LocalApproximation type " + type + " does not exist");
    }
}
