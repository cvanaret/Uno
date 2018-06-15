#include "GlobalizationMechanismFactory.hpp"
#include "TrustRegion.hpp"
#include "LineSearch.hpp"

std::shared_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(const std::string& type, GlobalizationStrategy& globalization_strategy) {
	if (type == "TR") {
		double radius = 10.;
		return std::make_shared<TrustRegion>(globalization_strategy, radius);
	}
	else if (type == "LS") {
		return std::make_shared<LineSearch>(globalization_strategy);
	}
	else {
		//TrustLineSearch globalization_mechanism = TrustLineSearch(*globalization_strategy, radius);
		throw std::invalid_argument("GlobalizationMechanism type " + type + " does not exist");
	}
}
