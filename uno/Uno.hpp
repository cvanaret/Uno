// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_H
#define UNO_H

#include "optimization/Result.hpp"
#include "optimization/TerminationStatus.hpp"

// forward declarations
class GlobalizationMechanism;
class Model;
class Options;
class Statistics;
class Timer;

/*
struct TimeOut : public std::exception {
   [[nodiscard]] const char* what() const noexcept override {
      return "The time limit was exceeded.\n";
   }
};
*/

class Uno {
public:
   Uno(GlobalizationMechanism& globalization_mechanism, const Options& options);

   [[nodiscard]] Result solve(const Model& model, Iterate& initial_iterate, const Options& options);

   static void print_uno_version();
   static void print_available_strategies();
   static void print_strategy_combination(const Options& options);
   static void print_optimization_summary(const Options& options, const Result& result);

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Globalization mechanism */
   const size_t max_iterations; /*!< Maximum number of iterations */
   const double time_limit; /*!< CPU time limit (can be inf) */

   void initialize(Statistics& statistics, Iterate& current_iterate, const Options& options);
   static Statistics create_statistics(const Model& model, const Options& options);
   [[nodiscard]] bool termination_criteria(TerminationStatus current_status, size_t iteration, double current_time) const;
   static void postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status);
   Result create_result(const Model& model, Iterate& current_iterate, size_t major_iterations, const Timer& timer);
};

#endif // UNO_H
