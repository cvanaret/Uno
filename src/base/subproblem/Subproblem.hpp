#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <cmath>
#include <vector>
#include <memory>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Phase.hpp"
#include "Direction.hpp"
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
    Subproblem(std::string residual_norm, std::vector<Range>& subproblem_variables_bounds, bool scale_residuals);
    virtual ~Subproblem();

    virtual Iterate evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
    
    virtual std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius=INFINITY) = 0;
    virtual std::vector<Direction> restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius=INFINITY) = 0;
    
    virtual void compute_optimality_measures(Problem& problem, Iterate& iterate) = 0;
    virtual void compute_infeasibility_measures(Problem& problem, Iterate& iterate, Direction& direction) = 0;
    
    static void project_point_in_bounds(std::vector<double>& x, const std::vector<Range>& variables_bounds);
    static double project_strictly_variable_in_bounds(double variable_value, const Range& variable_bounds);
    static std::vector<Range> generate_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints);
    static std::vector<double> compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, const std::vector<double>& default_multipliers, std::shared_ptr<LinearSolver> solver, double multipliers_max_size=1e3);
    static std::vector<double> compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, const std::vector<double>& default_multipliers, double multipliers_max_size=1e3);
    
    double compute_KKT_error(Problem& problem, Iterate& iterate, double objective_mutiplier);
    virtual double compute_complementarity_error_(Problem& problem, Iterate& iterate, Multipliers& multipliers);
    void compute_residuals(Problem& problem, Iterate& iterate, Multipliers& multipliers, double objective_multiplier);
    
    std::string residual_norm;
    // when the subproblem is reformulated (e.g. when slacks are introduced), the bounds may be altered as well
    std::vector<Range> subproblem_variables_bounds;
    int number_subproblems_solved;
    // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
    bool subproblem_definition_changed;
    bool scale_residuals;
};

#endif // SUBPROBLEM_H
