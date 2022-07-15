#include "GlobalizationStrategyFactory.hpp"
#include "MeritFunction.hpp"
#include "FilterStrategy.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
   if (strategy_type == "merit") {
      return std::make_unique<MeritFunction>(options);
   }
   else if (strategy_type == "filter" || strategy_type == "nonmonotone-filter") {
      return std::make_unique<FilterStrategy>(options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}
