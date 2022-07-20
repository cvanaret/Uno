// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

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

   virtual void initialize(Statistics& statistics, const Iterate& first_iterate) = 0;
   [[nodiscard]] virtual bool is_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
         double objective_multiplier, double predicted_reduction) = 0;
   virtual void reset() = 0;
   virtual void notify(Iterate& current_iterate) = 0;

   static void check_finiteness(const ProgressMeasures& progress);

protected:
   const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
   const double armijo_tolerance;

   [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
};

#endif // UNO_GLOBALIZATIONSTRATEGY_H
