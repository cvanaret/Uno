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

   const char* what() const throw() {
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
   InteriorPoint(Problem& problem, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region,
         bool scale_residuals);

   Iterate evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& default_multipliers) override;
   void evaluate_current_iterate(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) override;

   std::vector<Direction>
   compute_directions(Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   std::vector<Direction>
   restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) override;

   void compute_optimality_measures(const Problem& problem, Iterate& iterate) override;
   void compute_infeasibility_measures(const Problem& problem, Iterate& iterate, const Direction& direction) override;

   double constraint_violation(const Problem& problem, Iterate& iterate);
   double compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds);

   std::unique_ptr<HessianEvaluation> hessian_evaluation;
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that solves the subproblem */
   /* barrier parameter */
   double mu_optimality;
   double mu_feasibility;

   /* data structures */
   std::set<int> lower_bounded_variables; /* indices of the lower-bounded variables */
   std::set<int> upper_bounded_variables; /* indices of the upper-bounded variables */

private:
   void factorize_(COOMatrix& kkt_matrix, FunctionType problem_type);
   void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate);
   double evaluate_local_model_(Problem& problem, Iterate& current_iterate, std::vector<double>& solution);
   double barrier_function_(const Problem& problem, Iterate& iterate, const std::vector<Range>& variables_bounds);
   double compute_primal_length_(Iterate& iterate, std::vector<double>& ipm_solution, std::vector<Range>& variables_bounds, double tau);
   double
   compute_dual_length_(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z);
   COOMatrix assemble_optimality_kkt_matrix_(Problem& problem, Iterate& current_iterate);
   void modify_inertia_(COOMatrix& kkt_matrix, int size_first_block, int size_second_block, FunctionType problem_type);
   void generate_kkt_rhs_(Problem& problem, Iterate& current_iterate);
   std::vector<double> compute_lower_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution,
         std::vector<Range>& variables_bounds, double mu);
   std::vector<double> compute_upper_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution,
         std::vector<Range>& variables_bounds, double mu);
   Direction generate_direction_(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM);
   double compute_KKT_error_scaling_(Iterate& current_iterate) const;
   static double compute_predicted_reduction_(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length);

   bool force_symbolic_factorization = true;
   std::vector<double> rhs_;
   double inertia_hessian_;
   double inertia_hessian_last_;
   double inertia_constraints_;
   double default_multiplier_;
   size_t iteration_;
   size_t number_factorizations_;

   /* constants */
   InteriorPointParameters parameters_;
};

#endif // IPM_H
