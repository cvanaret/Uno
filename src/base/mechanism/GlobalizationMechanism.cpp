#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, int max_iterations)
      : globalization_strategy(globalization_strategy), max_iterations(max_iterations), number_iterations(0) {
}

std::optional<std::pair<Iterate, Direction> >
GlobalizationMechanism::find_first_acceptable_direction_(Statistics& statistics, Problem& problem, Iterate& current_iterate,
      std::vector<Direction>& directions, double step_length) {
   for (Direction& direction: directions) {
      try {
         std::optional<Iterate> acceptance_check =
               this->globalization_strategy.check_acceptance(statistics, problem, current_iterate, direction, step_length);
         if (acceptance_check.has_value()) {
            Iterate& accepted_iterate = acceptance_check.value();
            /* print summary */
            this->print_acceptance_();
            return std::pair<Iterate, Direction>(accepted_iterate, direction);
         }
      }
      catch (const NumericalError& e) {
      }
   }
   return std::nullopt;
}