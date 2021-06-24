#ifndef CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "ConstraintRelaxationStrategy.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static std::unique_ptr<ConstraintRelaxationStrategy> create(const std::string& constraint_relaxation_type, Problem& problem,
         const std::map<std::string, std::string>& options, bool use_trust_region);
};

#endif // CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
