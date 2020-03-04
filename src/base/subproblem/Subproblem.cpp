#include "Subproblem.hpp"

Subproblem::Subproblem() : number_subproblems_solved(0) {
}

Subproblem::~Subproblem() {
}

std::vector<Range> Subproblem::generate_constraints_bounds(Problem& problem, std::vector<double>& current_constraints) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb = problem.constraints_bounds[j].lb - current_constraints[j];
        double ub = problem.constraints_bounds[j].ub - current_constraints[j];
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}