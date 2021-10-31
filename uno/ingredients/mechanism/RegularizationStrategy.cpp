#include "RegularizationStrategy.hpp"

std::unique_ptr<RegularizationStrategy> RegularizationStrategyFactory::create() {
   return std::make_unique<InertiaBasedRegularization>();
}