#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(FeasibilityStrategy& feasibility_strategy, Subproblem& subproblem) :
   feasibility_strategy(feasibility_strategy), subproblem(subproblem) {
}