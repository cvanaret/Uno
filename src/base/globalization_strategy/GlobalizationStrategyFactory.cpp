#include "GlobalizationStrategyFactory.hpp"
#include "PenaltyStrategy.hpp"
#include "Filter.hpp"
#include "FilterStrategy.hpp"

std::shared_ptr<GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& type, LocalApproximation& local_approximation, double tolerance) {
	if (type == "penalty") {
		return std::make_shared<PenaltyStrategy>(local_approximation, tolerance);
	}
	else if (type == "filter") {
		TwoPhaseConstants constants = {0.1, 0.999};
		Tolerances filter_tolerances = {1e2, 1.25};
		
		/* create both filters */
		FilterConstants filter_constants = {0.999, 0.001};
		std::shared_ptr<Filter> filter_optimality = std::make_shared<Filter>(filter_constants);
		std::shared_ptr<Filter> filter_restoration = std::make_shared<Filter>(filter_constants);
		
		return std::make_shared<FilterStrategy>(local_approximation, filter_optimality, filter_restoration, constants, filter_tolerances, tolerance);
	}
	else if (type == "nonmonotone_filter") {
		TwoPhaseConstants constants = {0.1, 0.999};
		Tolerances filter_tolerances = {1e2, 1.25};
		
		/* create both filters */
		FilterConstants filter_constants = {0.999, 0.001};
		int number_dominated_entries = 3;
		std::shared_ptr<Filter> filter_optimality = std::make_shared<NonmonotoneFilter>(filter_constants, number_dominated_entries);
		std::shared_ptr<Filter> filter_restoration = std::make_shared<NonmonotoneFilter>(filter_constants, number_dominated_entries);
		
		return std::make_shared<FilterStrategy>(local_approximation, filter_optimality, filter_restoration, constants, filter_tolerances, tolerance);
	}
	else {
		throw std::invalid_argument("GlobalizationStrategy type " + type + " does not exist");
	}
}
