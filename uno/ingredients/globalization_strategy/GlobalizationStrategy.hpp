// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

#include "ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

/*! \class GlobalizationStrategy
 *  Ingredient that accepts or rejects a trial iterate
 */
class GlobalizationStrategy {
public:
   explicit GlobalizationStrategy(const Options& options);
   virtual ~GlobalizationStrategy() = default;

   virtual void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) = 0;
   [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) = 0;
   [[nodiscard]] virtual bool is_infeasibility_acceptable(const Model& model, Iterate& trial_iterate, Norm progress_norm) const = 0;

   virtual void reset() = 0;
   virtual void register_current_progress(const ProgressMeasures& current_progress) = 0;
   [[nodiscard]] virtual double get_infeasibility_upper_bound() const = 0;
   virtual void set_infeasibility_upper_bound(double new_upper_bound) = 0;

protected:
   const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
   const double armijo_tolerance;

   [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
};

#endif // UNO_GLOBALIZATIONSTRATEGY_H
