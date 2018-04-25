#include <iostream>
#include <cmath>
#include "AugmentedLagrangian.hpp"
#include "Utils.hpp"

AugmentedLagrangian::AugmentedLagrangian(QPSolver& solver): StepComputation("Augmented Lagrangian"), solver(solver) {
}

Step AugmentedLagrangian::compute_optimality_direction(Problem& problem, Iterate& current_point, double radius) {
	Step direction;

	return direction;
}

Step AugmentedLagrangian::compute_infeasibility_direction(Problem& problem, ConstraintPartition& constraint_partition, Iterate& current_point, double radius, std::vector<double>& d0) {
	Step direction = this->solver.solve(qp, d0);

	return direction;
}
