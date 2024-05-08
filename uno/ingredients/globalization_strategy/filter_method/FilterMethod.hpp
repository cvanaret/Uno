// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTERMETHOD_H
#define UNO_FILTERMETHOD_H

#include "../GlobalizationStrategy.hpp"
#include "filter/Filter.hpp"
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

class FilterMethod : public GlobalizationStrategy {
public:
   explicit FilterMethod(const Options& options);

   void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) override;
   void reset() override;
   void register_current_progress(const ProgressMeasures& current_progress) override;

protected:
   // pointer to allow polymorphism
   const std::unique_ptr<Filter> filter;
   const FilterStrategyParameters parameters; /*!< Set of constants */

   [[nodiscard]] static double unconstrained_merit_function(const ProgressMeasures& progress);
   [[nodiscard]] double compute_actual_objective_reduction(double current_objective_measure, double current_infeasibility, double trial_objective_measure);
   [[nodiscard]] bool switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const;
};

#endif // UNO_FILTERMETHOD_H
