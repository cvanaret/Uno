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
      const double upper_bound = stod(options.at("filter_ubd"));
      const double fact = stod(options.at("filter_fact"));
      FilterStrategyParameters strategy_parameters = {decrease_fraction, delta, upper_bound, fact};
      return std::make_unique<FilterStrategy>(strategy_parameters, options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}
