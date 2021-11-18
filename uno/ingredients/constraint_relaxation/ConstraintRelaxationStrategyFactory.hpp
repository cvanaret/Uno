#ifndef CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategyFactory {
public:
   static size_t get_max_number_variables(std::string_view constraint_relaxation_type, const Problem& problem);

   static std::unique_ptr<ConstraintRelaxationStrategy> create(const std::string& constraint_relaxation_type, Problem& problem, Subproblem&
      subproblem, const Options& options);
};

#endif // CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
