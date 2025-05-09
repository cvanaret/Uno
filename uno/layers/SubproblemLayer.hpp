// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMLAYER_H
#define UNO_SUBPROBLEMLAYER_H

#include <memory>

#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategyFactory.hpp"

namespace uno {
   // forward declarations
   class Options;
   class Statistics;

   class SubproblemLayer {
   public:
      std::unique_ptr<HessianModel> hessian_model;
      std::unique_ptr<RegularizationStrategy<double>> regularization_strategy;

      explicit SubproblemLayer(const Options& options):
         hessian_model(HessianModelFactory::create(options)),
         regularization_strategy(RegularizationStrategyFactory::create(options)) {
      }

      void initialize(const Model& model, const OptimizationProblem& problem) const {
         this->hessian_model->initialize(model);
         this->regularization_strategy->initialize_memory(problem, *this->hessian_model);
      }

      void initialize_statistics(Statistics& statistics, const Options& options) const {
         this->regularization_strategy->initialize_statistics(statistics, options);
      }

      [[nodiscard]] size_t get_hessian_evaluation_count() const {
         return this->hessian_model->evaluation_count;
      }
   };
} // namespace

#endif //UNO_SUBPROBLEMLAYER_H