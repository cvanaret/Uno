#ifndef IPM_H
#define IPM_H

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

        Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) override;

        SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) override;
        SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) override;
        SubproblemSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) override;

        void compute_measures(Problem& problem, Iterate& iterate) override;
        double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
        bool phase_1_required(SubproblemSolution& solution) override;
        
        double barrier_objective(Problem& problem, Iterate& iterate);
        double constraint_violation(Problem& problem, Iterate& iterate);
        double compute_central_complementarity_error(const Problem& problem, Iterate& iterate, double mu);
        
        MA57Solver solver; /*!< Solver that solves the subproblem */
        /* barrier parameter */
        double mu;

        /* data structures */
        std::vector<int> lower_bounded_variables; /* indices of the lower-bounded variables */
        std::vector<int> upper_bounded_variables; /* indices of the upper-bounded variables */
        std::map<int, int> lower_bounded_slacks; /* indices of the lower-bounded slacks */
        std::map<int, int> upper_bounded_slacks; /* indices of the upper-bounded slacks */
        MA57Factorization factorization;

    private:
        double evaluate_local_model(Problem& problem, Iterate& current_iterate, std::vector<double>& solution);
        std::map<int, double> barrier_function_gradient(Problem& problem, Iterate& current_iterate);
        double project_variable_in_bounds(double variable_value, Range& variable_bounds);
        std::vector<double> estimate_initial_multipliers(Problem& problem, Iterate& current_iterate);
        double compute_primal_length(Problem& problem, Iterate& iterate, std::vector<double>& ipm_solution, double tau, std::vector<Range>& variables_bounds);
        double compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z);
        COOMatrix generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds);
        void modify_inertia(COOMatrix& kkt_matrix, int number_variables, int number_constraints);
        std::vector<double> generate_kkt_rhs(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds);
        std::vector<double> compute_lower_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds);
        std::vector<double> compute_upper_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds);
        SubproblemSolution generate_direction(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM, std::vector<Range>& variables_bounds);
        double compute_KKT_error_scaling(Iterate& current_iterate);
        
        double inertia_hessian;
        double inertia_hessian_last;
        double inertia_constraints;
        
        /* constants */
        double tau_min;
        double default_multiplier;
        double k_sigma;
        double smax;
        double k_mu;
        double theta_mu;
        double k_epsilon;
        double multipliers_max_size;
        
        int iteration;
};

#endif // IPM_H
