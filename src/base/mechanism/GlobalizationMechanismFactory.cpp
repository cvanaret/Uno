#include <cmath>
#include "GlobalizationMechanismFactory.hpp"
#include "TrustRegion.hpp"
#include "LineSearch.hpp"
#include "TrustLineSearch.hpp"

std::shared_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(const std::string& type, GlobalizationStrategy& globalization_strategy, std::map<std::string, std::string> default_values) {
    if (type == "TR") {
        double radius = stod(default_values["TR_radius"]);
        int max_iterations = stoi(default_values["TR_max_iterations"]);
        return std::make_shared<TrustRegion>(globalization_strategy, radius, max_iterations);
    }
    else if (type == "LS") {
        int max_iterations = stoi(default_values["LS_max_iterations"]);
        double backtracking_ratio = stod(default_values["LS_backtracking_ratio"]);
        return std::make_shared<LineSearch>(globalization_strategy, max_iterations, backtracking_ratio);
    }
    else if (type == "TLS") {
        double radius = INFINITY;
        int max_iterations = stoi(default_values["LS_max_iterations"]);
        double ratio = stod(default_values["LS_ratio"]);
        return std::make_shared<TrustLineSearch>(globalization_strategy, radius, max_iterations, ratio);
    }
    else {
        throw std::invalid_argument("GlobalizationMechanism type " + type + " does not exist");
    }
}
