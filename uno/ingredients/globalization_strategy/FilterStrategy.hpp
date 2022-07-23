// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTERSTRATEGY_H
#define UNO_FILTERSTRATEGY_H

#include "GlobalizationStrategy.hpp"
#include "Filter.hpp"
#include "tools/Options.hpp"
#include "tools/Infinity.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for filter strategy
 *
 *  Set of constants to control the filter strategy
 */
struct FilterStrategyParameters {
   double delta; /*!< Switching constant */
   double upper_bound;
   double infeasibility_fraction;
   double switching_infeasibility_exponent;
};

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy : public GlobalizationStrategy {
public:
   explicit FilterStrategy(const Options& options);

   void initialize(const Iterate& first_iterate) override;
   bool is_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress, double objective_multiplier,
         double predicted_reduction) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   // pointer to allow polymorphism
   const std::unique_ptr<Filter> filter;
   double initial_filter_upper_bound{INF};
   const FilterStrategyParameters parameters; /*!< Set of constants */

   [[nodiscard]] bool switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const;
};

#endif // UNO_FILTERSTRATEGY_H
