#include "GlobalizationMechanismFactory.hpp"
#include "ingredients/mechanism/TrustRegion.hpp"
#include "ingredients/mechanism/BacktrackingLineSearch.hpp"

std::unique_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const
Options& options) {
   const std::string mechanism_type = options.at("mechanism");
    if (mechanism_type == "TR") {
        return std::make_unique<TrustRegion>(constraint_relaxation_strategy, options);
    }
    else if (mechanism_type == "LS") {
        return std::make_unique<BacktrackingLineSearch>(constraint_relaxation_strategy, options);
    }
    throw std::invalid_argument("GlobalizationMechanism " + mechanism_type + " is not supported");
}
