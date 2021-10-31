#ifndef REGULARIZATIONSTRATEGY_H
#define REGULARIZATIONSTRATEGY_H

#include <memory>

class RegularizationStrategy {
};

class InertiaBasedRegularization: public RegularizationStrategy {

};

class InertiaFreeRegularization: public RegularizationStrategy {

};

class RegularizationStrategyFactory {
public:
   static std::unique_ptr<RegularizationStrategy> create();
};

#endif // REGULARIZATIONSTRATEGY_H