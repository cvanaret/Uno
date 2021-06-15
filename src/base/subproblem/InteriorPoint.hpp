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

class InteriorPoint : public Subproblem {
public:
   InteriorPoint(const Problem& problem, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region,
         bool scale_residuals);

   Iterate evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& default_multipliers) override;
   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(const Problem& problem, Iterate& current_iterate) override;
   void compute_progress_measures(const Problem& problem, Iterate& iterate) override;
   int get_hessian_evaluation_count() const override;

private:
   /* barrier parameter */
   double barrier_parameter;
   std::unique_ptr<HessianEvaluation> hessian_evaluation;
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that solves the subproblem */
   /* constants */
   InteriorPointParameters parameters_;

   /* data structures */
   std::set<size_t> lower_bounded_variables; /* indices of the lower-bounded variables */
   std::set<size_t> upper_bounded_variables; /* indices of the upper-bounded variables */

   bool force_symbolic_factorization = true;
   std::vector<double> rhs_;
   double inertia_hessian_;
   double inertia_hessian_last_;
   double inertia_constraints_;
   double default_multiplier_;
   size_t iteration;
   size_t number_factorizations_;

   void factorize(COOMatrix& kkt_matrix, FunctionType problem_type);
   void evaluate_optimality_iterate(const Problem& problem, Iterate& current_iterate);
   static double compute_barrier_directional_derivative(const Problem&, Iterate& current_iterate, std::vector<double>& solution);
   double evaluate_barrier_function(const Problem& problem, Iterate& iterate, const std::vector<Range>& variables_bounds);
   double compute_primal_length(Iterate& current_iterate, std::vector<double>& ipm_solution, std::vector<Range>& variables_bounds, double tau);
   static double compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z);
   COOMatrix assemble_optimality_kkt_matrix(const Problem& problem, Iterate& current_iterate);
   void modify_inertia(COOMatrix& kkt_matrix, size_t size_first_block, size_t size_second_block, FunctionType problem_type);
   void generate_kkt_rhs(const Problem& problem, Iterate& current_iterate);
   std::vector<double> compute_lower_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution);
   std::vector<double> compute_upper_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution);
   Direction generate_direction(const Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM);
   double compute_KKT_error_scaling(Iterate& current_iterate) const;
   static double compute_predicted_reduction(Direction& direction, double step_length);
   static double constraint_violation(const Problem& problem, Iterate& iterate);
   double compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds) const;
};

#endif // IPM_H
