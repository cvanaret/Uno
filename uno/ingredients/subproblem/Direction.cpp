#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

Direction::Direction(size_t number_variables, size_t number_constraints):
   x(number_variables), multipliers(number_variables, number_constraints), constraint_partition(number_constraints) {

}

Direction::Direction(std::vector<double>& x, Multipliers& multipliers) : x(x), multipliers(multipliers), constraint_partition(multipliers
.constraints.size()) {
   const size_t number_variables = x.size();
   const size_t number_constraints = multipliers.constraints.size();
   this->active_set.bounds.at_lower_bound.reserve(number_variables);
   this->active_set.bounds.at_lower_bound.reserve(number_variables);
   this->active_set.constraints.at_lower_bound.reserve(number_constraints);
   this->active_set.constraints.at_upper_bound.reserve(number_constraints);
   this->constraint_partition.feasible.reserve(number_constraints);
   this->constraint_partition.infeasible.reserve(number_constraints);
}

std::string status_to_string(Status status) {
   switch (status) {
      case OPTIMAL:
         return "optimal";
      case UNBOUNDED_PROBLEM:
         return "unbounded problem";
      case BOUND_INCONSISTENCY:
         return "inconsistent bounds";
      case INFEASIBLE:
         return "infeasible";
      case INCORRECT_PARAMETER:
         return "incorrect parameter";
      case LP_INSUFFICIENT_SPACE:
         return "insufficient LP space";
      case HESSIAN_INSUFFICIENT_SPACE:
         return "insufficient Hessian space";
      case SPARSE_INSUFFICIENT_SPACE:
         return "insufficient sparse space";
      case MAX_RESTARTS_REACHED:
         return "max restarts reached";
      case UNDEFINED:
         return "undefined";
      default:
         return "unknown status, something went wrong";
   }
}

std::ostream& operator<<(std::ostream& stream, const Direction& direction) {
   stream << "Status: " << status_to_string(direction.status) << "\n";
   stream << "d^* = ";
   print_vector(stream, direction.x);

   stream << "evaluate_objective = " << direction.objective << "\n";
   stream << "norm = " << direction.norm << "\n";

   stream << "bound constraints active at lower bound =";
   for (size_t index: direction.active_set.bounds.at_lower_bound) {
      stream << " x" << index;
   }
   stream << "\n";
   stream << "bound constraints active at upper bound =";
   for (size_t index: direction.active_set.bounds.at_upper_bound) {
      stream << " x" << index;
   }
   stream << "\n";

   stream << "constraints at lower bound =";
   for (size_t index: direction.active_set.constraints.at_lower_bound) {
      stream << " c" << index;
   }
   stream << "\n";
   stream << "constraints at upper bound =";
   for (size_t index: direction.active_set.constraints.at_upper_bound) {
      stream << " c" << index;
   }
   stream << "\n";

   stream << "general feasible =";
   for (int j: direction.constraint_partition.feasible) {
      stream << " c" << j;
   }
   stream << "\n";

   stream << "general infeasible =";
   for (int j: direction.constraint_partition.infeasible) {
      stream << " c" << j;
      if (direction.constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         stream << " (lower)";
      }
      else if (direction.constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         stream << " (upper)";
      }
   }
   stream << "\n";

   stream << "lower bound multipliers = ";
   print_vector(stream, direction.multipliers.lower_bounds);
   stream << "upper bound multipliers = ";
   print_vector(stream, direction.multipliers.upper_bounds);
   stream << "constraint multipliers = ";
   print_vector(stream, direction.multipliers.constraints);

   return stream;
}
