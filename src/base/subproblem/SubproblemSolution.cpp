#include "SubproblemSolution.hpp"
#include "Logger.hpp"

SubproblemSolution::SubproblemSolution(std::vector<double>& x, Multipliers& multipliers):
x(x), multipliers(multipliers), objective_multiplier(1.), status(OPTIMAL), phase(OPTIMALITY), norm(0.), objective(0.),
is_descent_direction(true), constraint_partition(ConstraintPartition(multipliers.constraints.size())),
predicted_reduction([&](double step_length) { return 0.; }) {
}

std::ostream& operator<<(std::ostream &stream, SubproblemSolution& solution) {
    if (solution.status == OPTIMAL) {
        stream << "Status: optimal\n";
    }
    else if (solution.status == UNBOUNDED_PROBLEM) {
        stream << "Status: unbounded\n";
    }
    else if (solution.status == BOUND_INCONSISTENCY) {
        stream << "Status: bound inconsistency\n";
    }
    else if (solution.status == INFEASIBLE) {
        stream << "Status: infeasible subproblem\n";
    }
    else if (solution.status == INCORRECT_PARAMETER) {
        stream << "Status: incorrect parameter\n";
    }
    else if (solution.status == LP_INSUFFICIENT_SPACE) {
        stream << "Status: insufficient space for the LP\n";
    }
    else if (solution.status == HESSIAN_INSUFFICIENT_SPACE) {
        stream << "Status: insufficient space for the Hessian\n";
    }
    else if (solution.status == SPARSE_INSUFFICIENT_SPACE) {
        stream << "Status: insufficient space for the sparsity pattern\n";
    }
    else if (solution.status == MAX_RESTARTS_REACHED) {
        stream << "Status: maximum number of restarts reached\n";
    }
    else {
        stream << "Status " << solution.status << ": Beware peasant, something went wrong\n";
    }

    //stream << MAGENTA;
    stream << "d^* = ";
    print_vector(stream, solution.x);

    stream << "objective = " << solution.objective << "\n";
    stream << "norm = " << solution.norm << "\n";

    stream << "active set at upper bound =";
    for (unsigned int index: solution.active_set.at_upper_bound) {
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "active set at lower bound =";
    for (unsigned int index: solution.active_set.at_lower_bound) {
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "general feasible =";
    for (int j: solution.constraint_partition.feasible) {
        stream << " c" << j;
    }
    stream << "\n";

    stream << "general infeasible =";
    for (int j: solution.constraint_partition.infeasible) {
        stream << " c" << j;
        if (solution.constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            stream << " (lower)";
        }
        else if (solution.constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
            stream << " (upper)";
        }
    }
    stream << "\n";

    stream << "lower bound multipliers = ";
    print_vector(stream, solution.multipliers.lower_bounds);
    stream << "upper bound multipliers = ";
    print_vector(stream, solution.multipliers.upper_bounds);
    stream << "constraint multipliers = ";
    print_vector(stream, solution.multipliers.constraints);

    //stream << RESET;

    return stream;
}
