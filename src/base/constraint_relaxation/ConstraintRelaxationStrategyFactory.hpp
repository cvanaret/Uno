#ifndef CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "ConstraintRelaxationStrategy.hpp"
#include "Options.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static std::unique_ptr<ConstraintRelaxationStrategy> create(const std::string& constraint_relaxation_type, Problem& problem,
         const Options& options, bool use_trust_region);
};

#endif // CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
