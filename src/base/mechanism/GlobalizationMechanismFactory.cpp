#include <cmath>
#include "GlobalizationMechanismFactory.hpp"
#include "TrustRegion.hpp"
#include "LineSearch.hpp"
#include "TrustLineSearch.hpp"

std::shared_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(const std::string& type, GlobalizationStrategy& globalization_strategy, std::map<std::string, std::string> options) {
    double tolerance = stod(options["tolerance"]);

    if (type == "TR") {
        double radius = stod(options["TR_radius"]);
        int max_iterations = stoi(options["TR_max_iterations"]);
        return std::make_shared<TrustRegion>(globalization_strategy, tolerance, radius, max_iterations);
    }
    else if (type == "LS") {
        int max_iterations = stoi(options["LS_max_iterations"]);
        double backtracking_ratio = stod(options["LS_backtracking_ratio"]);
        return std::make_shared<LineSearch>(globalization_strategy, tolerance, max_iterations, backtracking_ratio);
    }
    else if (type == "TLS") {
        double radius = INFINITY;
        int max_iterations = stoi(options["LS_max_iterations"]);
        double ratio = stod(options["LS_ratio"]);
        return std::make_shared<TrustLineSearch>(globalization_strategy, tolerance, radius, max_iterations, ratio);
    }
    else {
        throw std::invalid_argument("GlobalizationMechanism type " + type + " does not exist");
    }
}
