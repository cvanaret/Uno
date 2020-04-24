#include "SubproblemSolution.hpp"
#include "Logger.hpp"

SubproblemSolution::SubproblemSolution(std::vector<double>& x, Multipliers& multipliers):
x(x), multipliers(multipliers), status(OPTIMAL), phase(OPTIMALITY), phase_1_required(false), norm(0.), objective(0.),
is_descent_direction(true), constraint_partition(ConstraintPartition(multipliers.constraints.size())) {
}

std::ostream& operator<<(std::ostream &stream, SubproblemSolution& solution) {
    unsigned int max_size = 50;

    if (solution.status == OPTIMAL) {
        //stream << GREEN "Status: optimal\n" RESET;
        stream << "Status: optimal\n";
    }
    else if (solution.status == UNBOUNDED_PROBLEM) {
        //stream << GREEN "Status: unbounded\n" RESET;
        stream << "Status: unbounded\n";
    }
    else {
        //stream << RED "Status " << solution.status << ": Beware peasant, something went wrong\n" RESET;
        stream << "Status " << solution.status << ": Beware peasant, something went wrong\n";
    }

    //stream << MAGENTA;
    stream << "d^* = ";
    print_vector(stream, solution.x, 0, max_size);

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
