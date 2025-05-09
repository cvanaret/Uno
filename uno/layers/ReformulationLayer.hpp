// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REFORMULATIONLAYER_H
#define UNO_REFORMULATIONLAYER_H

#include <memory>
#include <string>
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   // forward declarations
   class Options;
   class Statistics;

   class ReformulationLayer {
   public:
      std::unique_ptr<ConstraintRelaxationStrategy> constraint_relaxation_strategy;
      //std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      Direction direction{};

      ReformulationLayer(size_t number_constraints, size_t number_bounds_constraints, const Options& options):
         constraint_relaxation_strategy(ConstraintRelaxationStrategyFactory::create(number_constraints, number_bounds_constraints, options))
         // inequality_handling_method(InequalityHandlingMethodFactory::create(number_bounds_constraints, options))
         {
      }

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, const Options& options) {
         this->constraint_relaxation_strategy->initialize(statistics, model, initial_iterate, this->direction, options);
      }

      [[nodiscard]] std::string get_strategy_combination() const {
         return this->constraint_relaxation_strategy->get_name();
      }

      void compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
            Iterate& current_iterate, double trust_region_radius, WarmstartInformation& warmstart_information) {
         this->constraint_relaxation_strategy->compute_feasible_direction(statistics, globalization_strategy, model, current_iterate,
            this->direction, trust_region_radius, warmstart_information);
      }

      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
            Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
            WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const {
         return this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, globalization_strategy, model, current_iterate,
            trial_iterate, direction, step_length, warmstart_information, user_callbacks);
      }

      [[nodiscard]] size_t get_number_subproblems_solved() const {
         return this->constraint_relaxation_strategy->get_number_subproblems_solved();
      }
   };
} // namespace

#endif //UNO_REFORMULATIONLAYER_H