#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class BacktrackingLineSearch : public GlobalizationMechanism {
public:
   /*!
    *  Constructor
    */
   explicit BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations = 7, double backtracking_ratio
   = 0.5);

   Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
   std::tuple<Iterate, double, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;

   double step_length{1.};
   /* ratio of step length update in ]0, 1[ */
   const double backtracking_ratio;

private:
   const double min_step_length{1e-6};

   bool termination_();
   void print_iteration_();
   void add_statistics(Statistics& statistics, const Direction& direction);
   void update_step_length();

   /*
   double quadratic_interpolation_(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength = 1.);
   double cubic_interpolation_(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength1, double steplength2);
   double minimize_quadratic_(double a, double b);
   double minimize_cubic_(double a, double b, double c);
    */
};

#endif // LINESEARCH_H
