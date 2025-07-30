// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_H
#define UNO_H

#include <memory>
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Result.hpp"
#include "optimization/IterateStatus.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;
   class Statistics;
   class Timer;
   class UserCallbacks;

   class Uno {
   public:
      Uno(size_t number_constraints, const Options& options);

      // solve with or without user callbacks
      Result solve(const Model& model, Iterate& initial_iterate, const Options& options);
      Result solve(const Model& model, Iterate& initial_iterate, const Options& options, UserCallbacks& user_callbacks);

      static std::string current_version();
      static void print_available_strategies();
      void print_optimization_summary(const Result& result) const;

   private:
      std::unique_ptr<ConstraintRelaxationStrategy> constraint_relaxation_strategy;
      std::unique_ptr<GlobalizationStrategy> globalization_strategy;
      std::unique_ptr<GlobalizationMechanism> globalization_mechanism;
      Direction direction{};
      const size_t max_iterations; /*!< Maximum number of iterations */
      const double time_limit; /*!< CPU time limit (can be inf) */
      const bool print_solution;

      [[nodiscard]] std::string get_strategy_combination() const;
      void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, const Options& options);
      [[nodiscard]] static Statistics create_statistics(const Model& model, const Options& options);
      [[nodiscard]] bool termination_criteria(IterateStatus current_status, size_t iteration, double current_time,
         OptimizationStatus& optimization_status) const;
      static void postprocess_iterate(const Model& model, Iterate& iterate, IterateStatus termination_status);
      [[nodiscard]] Result create_result(const Model& model, OptimizationStatus optimization_status, Iterate& current_iterate,
         size_t major_iterations, const Timer& timer) const;
   };
} // namespace

#endif // UNO_H