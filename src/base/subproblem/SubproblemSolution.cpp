#include "SubproblemSolution.hpp"
#include "Logger.hpp"

//SubproblemSolution::SubproblemSolution(std::vector<double>& x, Multipliers& multipliers, Status& status, Phase& phase, bool phase_1_required, double norm, double objective, ConstraintActivity& active_set, ConstraintPartition& constraint_partition) :
//x(x), multipliers(multipliers), status(status), phase(phase), phase_1_required(phase_1_required), norm(norm), objective(objective), active_set(active_set), constraint_partition(constraint_partition) {
//}

SubproblemSolution::SubproblemSolution(std::vector<double>& x, Multipliers& multipliers, ActiveSet& active_set, ConstraintPartition& constraint_partition) :
x(x), multipliers(multipliers), phase_1_required(false), active_set(active_set), constraint_partition(constraint_partition) {
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
    for (unsigned int k = 0; k < solution.active_set.at_upper_bound.size(); k++) {
        unsigned int index = solution.active_set.at_upper_bound[k];
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "active set at lower bound =";
    for (unsigned int k = 0; k < solution.active_set.at_lower_bound.size(); k++) {
        unsigned int index = solution.active_set.at_lower_bound[k];
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "general feasible =";
    for (unsigned int k = 0; k < solution.constraint_partition.feasible_set.size(); k++) {
        int index = solution.constraint_partition.feasible_set[k];
        stream << " c" << index;
    }
    stream << "\n";

    stream << "general infeasible =";
    for (unsigned int k = 0; k < solution.constraint_partition.infeasible_set.size(); k++) {
        int index = solution.constraint_partition.infeasible_set[k];
        stream << " c" << index;
        if (solution.constraint_partition.constraint_status[index] == INFEASIBLE_LOWER) {
            stream << " (lower)";
        }
        else if (solution.constraint_partition.constraint_status[index] == INFEASIBLE_UPPER) {
            stream << " (upper)";
        }
    }
    stream << "\n";

    stream << "bound multipliers = ";
    print_vector(stream, solution.multipliers.bounds);

    stream << "constraint multipliers = ";
    print_vector(stream, solution.multipliers.constraints);

    //stream << RESET;

    return stream;
}
