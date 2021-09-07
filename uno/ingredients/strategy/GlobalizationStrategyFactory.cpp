#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
   if (strategy_type == "l1-penalty") {
      const double decrease_fraction = stod(options.at("decrease_fraction"));
      return std::make_unique<l1MeritFunction>(decrease_fraction);
   }
   else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
      const double decrease_fraction = stod(options.at("decrease_fraction"));
      const double delta = stod(options.at("filter_Delta"));
      const double ubd = stod(options.at("filter_ubd"));
      const double fact = stod(options.at("filter_fact"));
      FilterStrategyParameters strategy_parameters = {decrease_fraction, delta, ubd, fact};
      return std::make_unique<FilterStrategy>(strategy_parameters, options);
   }
   else {
      throw std::invalid_argument("GlobalizationStrategy type " + strategy_type + " does not exist");
   }
}
