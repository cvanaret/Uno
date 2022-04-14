#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include "ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static std::unique_ptr<ConstraintRelaxationStrategy> create(const Model& model, const Options& options);
};

#endif // UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
