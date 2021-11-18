#include "GlobalizationMechanismFactory.hpp"
#include "ingredients/mechanism/TrustRegion.hpp"
#include "ingredients/mechanism/BacktrackingLineSearch.hpp"

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
    throw std::invalid_argument("GlobalizationMechanism " + mechanism_type + " is not supported");
}
