#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"
#include "RegularizationStrategy.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class BacktrackingLineSearch : public GlobalizationMechanism {
public:
   explicit BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, const Scaling& scaling,
         Iterate& current_iterate) override;

private:
   std::unique_ptr<RegularizationStrategy> regularization_strategy;
   double step_length{1.};
   bool solving_feasibility_problem{false};
   /* ratio of step length update in ]0, 1[ */
   const double backtracking_ratio;
   const double min_step_length;
   const bool use_second_order_correction;

   Direction compute_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate);
   [[nodiscard]] bool termination() const;
   void print_iteration();
   void add_statistics(Statistics& statistics, const Direction& direction);
   void decrease_step_length();
};

#endif // UNO_BACKTRACKINGLINESEARCH_H
