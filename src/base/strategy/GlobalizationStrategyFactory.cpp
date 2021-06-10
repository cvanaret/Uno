#include "GlobalizationStrategyFactory.hpp"
#include "l1Penalty.hpp"
#include "Filter.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr<GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type,
      ConstraintRelaxationStrategy& feasibility_strategy, Subproblem& subproblem, const std::map<std::string, std::string>& options) {
    if (strategy_type == "penalty") {
        return std::make_unique<l1Penalty>(feasibility_strategy, subproblem);
    }
    else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
        double Sigma = stod(options.at("Sigma"));
        double Delta = stod(options.at("Delta"));
        double ubd = stod(options.at("ubd"));
        double fact = stod(options.at("fact"));
        FilterStrategyParameters strategy_parameters = {Sigma, Delta, ubd, fact};
        return std::make_unique<FilterStrategy>(feasibility_strategy, subproblem, strategy_parameters, options);
    }
    else {
        throw std::invalid_argument("GlobalizationStrategy type " + strategy_type + " does not exist");
    }
}
