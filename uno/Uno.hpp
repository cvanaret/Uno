// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_H
#define UNO_H

#include "optimization/Model.hpp"
#include "optimization/Result.hpp"
#include "optimization/TerminationStatus.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanism.hpp"

class Uno {
public:
   Uno(GlobalizationMechanism& globalization_mechanism, const Options& options);

   [[nodiscard]] Result solve(Statistics& statistics, const Model& model, Iterate& initial_iterate);
   static void print_available_strategies();

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
   const double tolerance; /*!< Tolerance of the termination criteria */
   const size_t max_iterations; /*!< Maximum number of iterations */
   const bool terminate_with_small_step;
   const double small_step_threshold;

   static void add_statistics(Statistics& statistics, const Model& model, const Iterate& iterate, size_t major_iterations);
   [[nodiscard]] bool termination_criterion(TerminationStatus current_status, size_t iteration) const;
   [[nodiscard]] TerminationStatus check_termination(const Model& model, Iterate& current_iterate, double step_norm) const;
   static void postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status);
};

#endif // UNO_H
