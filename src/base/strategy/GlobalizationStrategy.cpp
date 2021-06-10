#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, Subproblem& subproblem) :
   constraint_relaxation_strategy(constraint_relaxation_strategy), subproblem(subproblem) {
}