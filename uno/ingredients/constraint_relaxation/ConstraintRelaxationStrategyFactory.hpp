#ifndef CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static std::unique_ptr<ConstraintRelaxationStrategy> create(const std::string& constraint_relaxation_type, Problem& problem,
         const Options& options);
};

#endif // CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
