#include "SubproblemFactory.hpp"
#include "QPApproximation.hpp"
#include "BQPDSolver.hpp"
#include "InteriorPoint.hpp"

std::shared_ptr<Subproblem> SubproblemFactory::create(const std::string& type, QPSolver& solver, std::map<std::string, std::string> default_values) {
    if (type == "QP") {
        return std::make_shared<QPApproximation>(solver);
    }
    else if (type == "IPM") {
        return std::make_shared<InteriorPoint>();
    }
    else {
        throw std::invalid_argument("LocalApproximation type " + type + " does not exist");
    }
}
