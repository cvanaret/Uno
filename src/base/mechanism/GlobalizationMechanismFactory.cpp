#include "GlobalizationMechanismFactory.hpp"
#include "TrustRegion.hpp"
#include "BacktrackingLineSearch.hpp"
//#include "DualUpdate.hpp"
//#include "TrustLineSearch.hpp"

std::unique_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(const std::string& mechanism_type, ConstraintRelaxationStrategy&
constraint_relaxation_strategy, const Options& options) {
    if (mechanism_type == "TR") {
        double radius = stod(options.at("TR_radius"));
        int max_iterations = std::stoi(options.at("TR_max_iterations"));
        return std::make_unique<TrustRegion>(constraint_relaxation_strategy, radius, max_iterations);
    }
    else if (mechanism_type == "LS") {
        int max_iterations = std::stoi(options.at("LS_max_iterations"));
        double backtracking_ratio = std::stod(options.at("LS_backtracking_ratio"));
        return std::make_unique<BacktrackingLineSearch>(constraint_relaxation_strategy, max_iterations, backtracking_ratio);
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
        throw std::invalid_argument("GlobalizationMechanism type " + mechanism_type + " does not exist");
    }
}
