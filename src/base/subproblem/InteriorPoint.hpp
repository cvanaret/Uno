#ifndef IPM_H
#define IPM_H

#include <exception>
#include <set>
#include "Subproblem.hpp"
#include "LinearSolver.hpp"
#include "HessianEvaluation.hpp"

struct InteriorPointParameters {
    double tau_min;
    double k_sigma;
    double smax;
    double k_mu;
    double theta_mu;
    double k_epsilon;
    double kappa;
};

struct UnstableInertiaCorrection : public std::exception {

    const char* what() const throw () {
        return "The inertia correction got unstable (delta_w > 1e40)";
    }
};

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
    InteriorPoint(Problem& problem, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);

    Iterate evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& default_multipliers) override;

    SubproblemSolution compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY) override;
    SubproblemSolution restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY) override;

    void compute_optimality_measures(Problem& problem, Iterate& iterate) override;
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) override;
    
    double constraint_violation(Problem& problem, Iterate& iterate);
    double compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds);

    std::shared_ptr<HessianEvaluation> hessian_evaluation;
    std::shared_ptr<LinearSolver> solver; /*!< Solver that solves the subproblem */
    /* barrier parameter */
    double mu_optimality;
    double mu_feasibility;

    /* data structures */
    std::set<int> lower_bounded_variables; /* indices of the lower-bounded variables */
    std::set<int> upper_bounded_variables; /* indices of the upper-bounded variables */

private:
    void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate);
    double evaluate_local_model_(Problem& problem, Iterate& current_iterate, std::vector<double>& solution);
    double barrier_function_(Problem& problem, Iterate& iterate, std::vector<Range>& variables_bounds);
    double project_variable_in_bounds_(double variable_value, Range& variable_bounds);
    double compute_primal_length_(Iterate& iterate, std::vector<double>& ipm_solution, std::vector<Range>& variables_bounds, double tau);
    double compute_dual_length_(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z);
    COOMatrix generate_optimality_kkt_matrix_(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds);
    void modify_inertia_(COOMatrix& kkt_matrix, int size_first_block, int size_second_block);
    std::vector<double> generate_kkt_rhs_(Problem& problem, Iterate& current_iterate);
    std::vector<double> compute_lower_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds, double mu);
    std::vector<double> compute_upper_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds, double mu);
    SubproblemSolution generate_direction_(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM);
    double compute_KKT_error_scaling_(Iterate& current_iterate);
    double compute_predicted_reduction_(SubproblemSolution& solution, double step_length);

    double inertia_hessian_;
    double inertia_hessian_last_;
    double inertia_constraints_;
    double default_multiplier_;
    int iteration_;

    /* constants */
    InteriorPointParameters parameters_;
};

#endif // IPM_H
