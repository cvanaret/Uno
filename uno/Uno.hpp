// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_H
#define UNO_H

#include "optimization/Model.hpp"
#include "optimization/Result.hpp"
#include "optimization/TerminationStatus.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanism.hpp"

struct TimeOut : public std::exception {
   [[nodiscard]] const char* what() const noexcept override {
      return "The time limit was exceeded.\n";
   }
};

class Uno {
public:
   Uno(GlobalizationMechanism& globalization_mechanism, const Options& options);

   [[nodiscard]] Result solve(Statistics& statistics, const Model& model, Iterate& initial_iterate);
   static void print_available_strategies();

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Globalization mechanism */
   const size_t max_iterations; /*!< Maximum number of iterations */
   const double time_limit; /*!< CPU time limit (can be inf) */

   static void add_statistics(Statistics& statistics, const Iterate& iterate, size_t major_iterations);
   [[nodiscard]] bool termination_criteria(TerminationStatus current_status, size_t iteration, double current_time) const;
   static void postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status);
};

#endif // UNO_H
