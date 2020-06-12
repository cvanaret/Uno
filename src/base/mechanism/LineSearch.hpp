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
    LineSearch(GlobalizationStrategy& globalization_strategy, double tolerance, int max_iterations = 30, double backtracking_ratio = 0.5);

    Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
    Iterate compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) override;

    double step_length;
    /* ratio of step length update in ]0, 1[ */
    double backtracking_ratio;

private:
    double min_step_length;

    std::vector<Range> compute_subproblem_bounds_(Iterate& current_iterate);
    bool termination_(bool is_accepted);
    void print_iteration_();
    void print_acceptance_() override;
    void print_warning_(const char* message);

    double quadratic_interpolation_(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength = 1.);
    double cubic_interpolation_(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength1, double steplength2);
    double minimize_quadratic_(double a, double b);
    double minimize_cubic_(double a, double b, double c);
};

#endif // LINESEARCH_H
