#include "GlobalizationStrategyFactory.hpp"
#include "PenaltyMeritFunction.hpp"
#include "Filter.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr<GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, Subproblem& subproblem, std::map<std::string, std::string>& options) {
    if (strategy_type == "penalty") {
        return std::make_unique<PenaltyMeritFunction>(subproblem);
    }
    else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
        double Sigma = stod(options["Sigma"]);
        double Delta = stod(options["Delta"]);
        double ubd = stod(options["ubd"]);
        double fact = stod(options["fact"]);
        FilterStrategyParameters strategy_parameters = {Sigma, Delta, ubd, fact};
        return std::make_unique<FilterStrategy>(subproblem, strategy_parameters, options);
    }
    else {
        throw std::invalid_argument("GlobalizationStrategy type " + strategy_type + " does not exist");
    }
}
