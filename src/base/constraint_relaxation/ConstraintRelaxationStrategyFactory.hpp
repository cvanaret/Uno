#ifndef CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "ConstraintRelaxationStrategy.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static std::unique_ptr<ConstraintRelaxationStrategy> create(const std::string& constraint_relaxation_type, Problem& problem, Subproblem&
   subproblem, const std::map<std::string, std::string>& options);
};

#endif // CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
