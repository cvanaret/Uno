#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(Subproblem& subproblem, double tolerance) : subproblem(subproblem), tolerance(tolerance) {
}

GlobalizationStrategy::~GlobalizationStrategy() {
}