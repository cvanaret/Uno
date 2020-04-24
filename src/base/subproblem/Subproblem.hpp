#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <cmath>
#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Phase.hpp"
#include "SubproblemSolution.hpp"
#include "Constraint.hpp"
#include "MA57Solver.hpp"

/*! \class Subproblem
 * \brief Subproblem
 *
 *  Local approximation of a nonlinear optimization problem (virtual class) 
 */
class Subproblem {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     * \param name: name of the strategy
     */
    Subproblem(std::string residual_norm);
    virtual ~Subproblem();

    virtual Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, bool use_trust_region) = 0;
    
    // TODO return a list of steps
    virtual SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius=INFINITY) = 0;
    virtual SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius=INFINITY) = 0;
    
    virtual void compute_optimality_measures(Problem& problem, Iterate& iterate) = 0;
    virtual void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) = 0;
    
    virtual double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) = 0;
    virtual bool phase_1_required(SubproblemSolution& solution) = 0;

    std::vector<Range> generate_variables_bounds(Iterate& current_iterate, double trust_region_radius);
    static double project_variable_in_bounds(double variable_value, Range& variable_bounds);
    static std::vector<Range> generate_constraints_bounds(Problem& problem, std::vector<double>& current_constraints);
    static std::vector<double> compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, std::vector<double>& default_multipliers, MA57Solver& solver, double multipliers_max_size=1e3);
    static std::vector<double> compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, std::vector<double>& default_multipliers, double multipliers_max_size=1e3);
    
    std::string residual_norm;
    // when the subproblem is reformulated (e.g. when slacks are introduced), the bounds may be altered as well
    std::vector<Range> subproblem_variables_bounds;
    int number_subproblems_solved;
    // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
    bool subproblem_definition_changed;
};

#endif // SUBPROBLEM_H
