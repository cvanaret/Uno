#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
#include "Filter.hpp"
#include "tools/Options.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for filter strategy
 *
 *  Set of constants to control the filter strategy
 */
struct FilterStrategyParameters {
   double decrease_fraction; /*!< Sufficient reduction constant */
   double Delta; /*!< Switching constant */
   double upper_bound;
   double infeasibility_factor;
};

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy : public GlobalizationStrategy {
public:
   FilterStrategy(FilterStrategyParameters strategy_constants, const Options& options);

   void initialize(Statistics& statistics, const Iterate& first_iterate) override;
   bool check_acceptance(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
         double objective_multiplier, double predicted_reduction) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   // use pointers to allow polymorphism
   const std::unique_ptr<Filter> filter;
   double initial_filter_upper_bound;
   const FilterStrategyParameters parameters; /*!< Set of constants */

   static bool switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction);
   static bool armijo_condition(double predicted_reduction, double actual_reduction, double decrease_fraction);
};

#endif // FILTERSTRATEGY_H
