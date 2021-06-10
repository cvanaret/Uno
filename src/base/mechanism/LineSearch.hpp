#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class LineSearch : public GlobalizationMechanism {
public:
   /*!
    *  Constructor
    */
   explicit LineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations = 30, double backtracking_ratio = 0.5);

   Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
   std::pair<Iterate, Direction> compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) override;

   double step_length;
   /* ratio of step length update in ]0, 1[ */
   double backtracking_ratio;

private:
   double min_step_length;

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
