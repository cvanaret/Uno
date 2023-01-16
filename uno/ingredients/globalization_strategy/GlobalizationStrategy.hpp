// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

#include "PredictedReductionModel.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

/*! \class GlobalizationStrategy
 * \brief Step acceptance strategy
 *
 *  Strategy that accepts or declines a trial step (virtual class)
 */
class GlobalizationStrategy {
public:
   explicit GlobalizationStrategy(const Options& options);
   virtual ~GlobalizationStrategy() = default;

   virtual void initialize(const Iterate& first_iterate) = 0;
   [[nodiscard]] virtual bool is_iterate_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
         const PredictedReduction& predicted_reduction, double objective_multiplier) = 0;
   [[nodiscard]] virtual bool is_infeasibility_acceptable(double infeasibility_measure) const = 0;

   virtual void reset() = 0;
   virtual void register_current_progress(const ProgressMeasures& current_progress) = 0;

protected:
   const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
   const double armijo_tolerance;

   [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
   static void check_finiteness(const ProgressMeasures& progress, double objective_multiplier);
};

#endif // UNO_GLOBALIZATIONSTRATEGY_H
