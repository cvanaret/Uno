#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr<GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
    if (strategy_type == "l1-penalty") {
        return std::make_unique<l1MeritFunction>();
    }
    else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
        double Sigma = stod(options.at("Sigma"));
        double Delta = stod(options.at("Delta"));
        double ubd = stod(options.at("ubd"));
        double fact = stod(options.at("fact"));
        FilterStrategyParameters strategy_parameters = {Sigma, Delta, ubd, fact};
        return std::make_unique<FilterStrategy>(strategy_parameters, options);
    }
    else {
        throw std::invalid_argument("GlobalizationStrategy type " + strategy_type + " does not exist");
    }
}
