#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(Subproblem& subproblem, double tolerance) : subproblem(subproblem) {
    this->tolerance = tolerance;
}

GlobalizationStrategy::~GlobalizationStrategy() {
}