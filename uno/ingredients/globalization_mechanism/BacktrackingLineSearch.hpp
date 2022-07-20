// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class BacktrackingLineSearch : public GlobalizationMechanism {
public:
   BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Statistics& statistics, Iterate& first_iterate) override;
   [[nodiscard]] std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, Iterate& current_iterate) override;

private:
   double step_length{1.};
   bool solving_feasibility_problem{false};
   // ratio of step length update in ]0, 1[
   const double backtracking_ratio;
   const double min_step_length;
   const bool use_second_order_correction;
   // statistics table
   const int statistics_SOC_column_order;
   const int statistics_LS_step_length_column_order;

   [[nodiscard]] Direction compute_direction(Statistics& statistics, Iterate& current_iterate);
   [[nodiscard]] bool termination() const;
   void print_iteration();
   void add_statistics(Statistics& statistics, const Direction& direction);
   void decrease_step_length();
};

#endif // UNO_BACKTRACKINGLINESEARCH_H
