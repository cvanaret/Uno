#ifndef INTERIORPOINT_H
#define INTERIORPOINT_H

#include "Subproblem.hpp"
#include "MA57Solver.hpp"

/*! \class InteriorPoint
 * \brief Interior Point Method
 *
 *  Implementation of an Interior Point Method
 */
class InteriorPoint : public Subproblem {
    public:
        /*!
         *  Constructor
         */
        InteriorPoint();

        Iterate initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, int number_variables, int number_constraints, bool use_trust_region) override;

        LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) override;
        LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) override;
        LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) override;

        void compute_measures(Problem& problem, Iterate& iterate) override;
        double objective(Problem& problem, Iterate& iterate);
        double constraint_violation(Problem& problem, Iterate& iterate);
        double evaluate_local_model(Problem& problem, Iterate& current_iterate, std::vector<double>& solution);
        
        MA57Solver solver; /*!< Solver that solves the subproblem */

        /* barrier parameter */
        double mu;

        /* data structures */
        std::vector<int> lower_bounded_variables; /* indices of the variables with lower bounds */
        std::vector<int> upper_bounded_variables; /* indices of the variables with upper bounds */
        std::vector<int> slacked_constraints; /* indices of the inequality constraints that need a slack variable */
        std::vector<int> lower_bounded_slacks; /* indices of the slacks with lower bounds */
        std::vector<int> upper_bounded_slacks; /* indices of the slacks with upper bounds */
        MA57Data data;

    private:
        double project_variable_in_bounds(double current_value, double lb, double ub);
        std::vector<double> estimate_initial_multipliers(Problem& problem, Iterate& current_iterate);
        double compute_primal_length(Problem& problem, Iterate& iterate, std::vector<double>& ipm_solution, double tau, std::vector<double> variable_lb, std::vector<double> variable_ub);
        double compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& delta_z);
        COOMatrix generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
        void inertia_correction(Problem& problem, COOMatrix& kkt_matrix);
        std::vector<double> generate_rhs(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
        std::vector<double> compute_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<double>& variable_lb, std::vector<double>& variable_ub);
        double update_barrier_parameter(Problem& problem, Iterate& current_iterate);
        double compute_error(Problem& problem, Iterate& iterate);
        
        double inertia_hessian;
        double inertia_hessian_last;
        double inertia_constraints;
        
        /* constants */
        double tau_min;
        double default_multiplier;
        double k_sigma;
        
        int iteration;
};

#endif // INTERIORPOINT_H
