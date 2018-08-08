#include "GlobalizationStrategyFactory.hpp"
#include "PenaltyStrategy.hpp"
#include "Filter.hpp"
#include "FilterStrategy.hpp"

std::shared_ptr<GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& type, Subproblem& subproblem,
		std::map<std::string, std::string> default_values) {
	double tolerance = stod(default_values["tolerance"]);
	
	if (type == "penalty") {
		return std::make_shared<PenaltyStrategy>(subproblem, tolerance);
	}
	else if (type == "filter" || type == "nonmonotone_filter") {
		double Sigma = stod(default_values["Sigma"]);
		double Delta = stod(default_values["Delta"]);
		double ubd = stod(default_values["ubd"]);
		double fact = stod(default_values["fact"]);
		TwoPhaseConstants constants = {Sigma, Delta, ubd, fact};
		
		/* create both filters */
		double Beta = stod(default_values["Beta"]);
		double Gamma = stod(default_values["Gamma"]);
		FilterConstants filter_constants = {Beta, Gamma};
		if (type == "filter") {
			std::shared_ptr<Filter> filter_optimality = std::make_shared<Filter>(filter_constants);
			std::shared_ptr<Filter> filter_restoration = std::make_shared<Filter>(filter_constants);
			return std::make_shared<FilterStrategy>(subproblem, filter_optimality, filter_restoration, constants, tolerance);
		}
		else {
			int number_dominated_entries = stoi(default_values["number_dominated_entries"]);
			std::shared_ptr<Filter> filter_optimality = std::make_shared<NonmonotoneFilter>(filter_constants, number_dominated_entries);
			std::shared_ptr<Filter> filter_restoration = std::make_shared<NonmonotoneFilter>(filter_constants, number_dominated_entries);
			return std::make_shared<FilterStrategy>(subproblem, filter_optimality, filter_restoration, constants, tolerance);
		}
	}
	else {
		throw std::invalid_argument("GlobalizationStrategy type " + type + " does not exist");
	}
}
