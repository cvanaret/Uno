#include "LocalSolution.hpp"
#include "Logger.hpp"

LocalSolution::LocalSolution(std::vector<double>& x, int n, int m):
		x(x), multipliers(n + m) {
}

std::ostream& operator<< (std::ostream &stream, LocalSolution& solution) {
	unsigned int max_size = 50;
	
	if (solution.status == OPTIMAL) {
		stream << GREEN "Status: optimal\n" RESET;
	}
	else if (solution.status == UNBOUNDED) {
		stream << GREEN "Status: unbounded\n" RESET;
	}
	else {
		stream << RED "Status " << solution.status << ": Beware peasant, something went wrong: " << solution.status << "\n" RESET;
	}
	
	stream << MAGENTA;
	stream << "x^* = ";
	print_vector(stream, solution.x, max_size);
	
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
	}
	stream << "\n";

	stream << "multipliers = ";
	print_vector(stream, solution.multipliers);
	
	stream << RESET;

	return stream;
}
