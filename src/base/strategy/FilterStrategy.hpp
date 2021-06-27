#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
#include "Filter.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for filter and tube strategies
 *
 *  Set of constants to control the filter and tube strategies
 */
struct FilterStrategyParameters {
   double Sigma; /*!< Sufficient reduction constant */
   double Delta; /*!< Switching constant */
   double ubd;
   double fact;
};

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy : public GlobalizationStrategy {
public:
   FilterStrategy(FilterStrategyParameters& strategy_constants, const std::map<std::string, std::string>& options);

   /* use pointers to allow polymorphism */
   const std::unique_ptr<Filter> filter;
   double initial_filter_upper_bound;

   void initialize(Statistics& statistics, const Iterate& first_iterate) override;
   bool check_acceptance(Statistics& statistics, ProgressMeasures& current_progress, ProgressMeasures& trial_progress, double objective_multiplier,
         double predicted_reduction) override;
   void reset() override;
   void notify(Iterate& current_iterate) override;

private:
   const FilterStrategyParameters parameters_; /*!< Set of constants */
};

#endif // FILTERSTRATEGY_H
