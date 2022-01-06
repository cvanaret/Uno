#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
   if (strategy_type == "l1-merit") {
      return std::make_unique<l1MeritFunction>(options);
   }
   else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
      return std::make_unique<FilterStrategy>(options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}
