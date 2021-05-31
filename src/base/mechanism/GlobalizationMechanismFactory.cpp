#include "GlobalizationMechanismFactory.hpp"
#include "TrustRegion.hpp"
#include "LineSearch.hpp"
//#include "DualUpdate.hpp"
//#include "TrustLineSearch.hpp"

std::unique_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(const std::string& type, GlobalizationStrategy&
globalization_strategy, std::map<std::string, std::string>& options) {
    double tolerance = std::stod(options["tolerance"]);

    if (type == "TR") {
        double radius = stod(options["TR_radius"]);
        int max_iterations = std::stoi(options["TR_max_iterations"]);
        return std::make_unique<TrustRegion>(globalization_strategy, tolerance, radius, max_iterations);
    }
    else if (type == "LS") {
        int max_iterations = std::stoi(options["LS_max_iterations"]);
        double backtracking_ratio = std::stod(options["LS_backtracking_ratio"]);
        return std::make_unique<LineSearch>(globalization_strategy, tolerance, max_iterations, backtracking_ratio);
    }
//    else if (type == "TLS") {
//        double radius = INFINITY;
//        int max_iterations = std::stoi(options["LS_max_iterations"]);
//        double ratio = std::stod(options["LS_ratio"]);
//        return std::make_unique<TrustLineSearch>(globalization_strategy, tolerance, radius, max_iterations, ratio);
//    }
    //else if (type == "dual_update") {
    //    return std::make_unique<DualUpdate>(globalization_strategy, tolerance);
    //}
    else {
        throw std::invalid_argument("GlobalizationMechanism type " + type + " does not exist");
    }
}
