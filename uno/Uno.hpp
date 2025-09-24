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
#include "optimization/SolutionStatus.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;
   class Statistics;
   class Timer;
   class UserCallbacks;

   class Uno {
   public:
      Uno() = default;

      // solve with or without user callbacks
      Result solve(const Model& model, const Options& options);
      Result solve(const Model& model, const Options& options, UserCallbacks& user_callbacks);

      static std::string current_version();
      static void print_available_strategies();

   private:
      std::unique_ptr<ConstraintRelaxationStrategy> constraint_relaxation_strategy{};
      std::unique_ptr<GlobalizationStrategy> globalization_strategy{};
      std::unique_ptr<GlobalizationMechanism> globalization_mechanism{};
      Direction direction{};

      void pick_ingredients(const Model& model, const Options& options);
      void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, const Options& options);
      [[nodiscard]] static Statistics create_statistics(const Model& model, const Options& options);
      [[nodiscard]] static bool termination_criteria(SolutionStatus solution_status, size_t iteration, size_t max_iterations,
         double current_time, double time_limit, OptimizationStatus& optimization_status);
      [[nodiscard]] Result uno_solve(const Model& model, const Options& options, UserCallbacks& user_callbacks);
      static void postprocess_iterate(const Model& model, Iterate& iterate);
      [[nodiscard]] Result create_result(const Model& model, OptimizationStatus optimization_status, Iterate& solution,
         size_t major_iterations, const Timer& timer) const;
      [[nodiscard]] std::string get_strategy_combination() const;
      void print_optimization_summary(const Result& result, bool print_solution) const;
   };
} // namespace

#endif // UNO_H