// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREDICTEDREDUCTIONMODELS_H
#define UNO_PREDICTEDREDUCTIONMODELS_H

#include <functional>
#include <utility>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"

namespace uno {
   // lazy local models of the predicted reductions. For a fixed (current iterate, direction), the
   // step-independent quantities (Jacobian-direction product, directional derivatives, dᵀHd, ...) are
   // evaluated once when the model is built; each closure below only assembles the reduction for a given
   // step length (and, for the objective, objective parameter)
   using PredictedInfeasibilityReduction = std::function<double(double step_length)>;
   using PredictedObjectiveReduction = std::function<double(double step_length, double objective_multiplier,
      bool use_curvature_information)>;
   using PredictedAuxiliaryReduction = std::function<double(double step_length)>;

   // built once per (subproblem, current iterate, direction); queried once by the trust-region method, and repeatedly
   // by the backtracking line search
   class PredictedReductionModels {
   public:
      PredictedReductionModels(PredictedInfeasibilityReduction infeasibility,
            PredictedObjectiveReduction objective, PredictedAuxiliaryReduction auxiliary):
         infeasibility_model(std::move(infeasibility)),
         objective_model(std::move(objective)),
         auxiliary_model(std::move(auxiliary)) {
      }

      // assemble the predicted reductions for a given step length
      [[nodiscard]] ProgressMeasures operator()(double step_length) const {
         return {
            this->infeasibility_model(step_length),
            // predicted objective reduction left as a function of the objective multiplier
            [objective = this->objective_model, step_length](double objective_multiplier) {
               return objective(step_length, objective_multiplier, /* use_curvature_information = */ true);
            },
            this->auxiliary_model(step_length)
         };
      }

   private:
      PredictedInfeasibilityReduction infeasibility_model;
      PredictedObjectiveReduction objective_model;
      PredictedAuxiliaryReduction auxiliary_model;
   };
} // namespace

#endif // UNO_PREDICTEDREDUCTIONMODELS_H