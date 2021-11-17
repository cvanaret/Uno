#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

Direction::Direction(size_t number_variables, size_t number_constraints):
   x(number_variables), multipliers(number_variables, number_constraints) {
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

   stream << "objective = " << direction.objective << "\n";
   stream << "norm = " << direction.norm << "\n";

   stream << "bound constraints active at lower bound =";
   for (size_t i: direction.active_set.bounds.at_lower_bound) {
      stream << " x" << i;
   }
   stream << "\n";
   stream << "bound constraints active at upper bound =";
   for (size_t i: direction.active_set.bounds.at_upper_bound) {
      stream << " x" << i;
   }
   stream << "\n";

   stream << "constraints at lower bound =";
   for (size_t j: direction.active_set.constraints.at_lower_bound) {
      stream << " c" << j;
   }
   stream << "\n";
   stream << "constraints at upper bound =";
   for (size_t j: direction.active_set.constraints.at_upper_bound) {
      stream << " c" << j;
   }
   stream << "\n";

   if (direction.constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      stream << "general feasible =";
      for (size_t j: constraint_partition.feasible) {
         stream << " c" << j;
      }
      stream << "\n";

      stream << "\ngeneral lower infeasible =";
      for (size_t j: constraint_partition.lower_bound_infeasible) {
         stream << " c" << j << " ";
      }
      stream << "\ngeneral upper infeasible =";
      for (size_t j: constraint_partition.upper_bound_infeasible) {
         stream << " c" << j << " ";
      }
      stream << "\n";
   }

   stream << "lower bound multipliers = ";
   print_vector(stream, direction.multipliers.lower_bounds);
   stream << "upper bound multipliers = ";
   print_vector(stream, direction.multipliers.upper_bounds);
   stream << "constraint multipliers = ";
   print_vector(stream, direction.multipliers.constraints);

   return stream;
}
