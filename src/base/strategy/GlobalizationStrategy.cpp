#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(ConstraintRelaxationStrategy& feasibility_strategy, Subproblem& subproblem) :
   feasibility_strategy(feasibility_strategy), subproblem(subproblem) {
}