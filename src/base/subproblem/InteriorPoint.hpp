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
   InteriorPoint(const Problem& problem, const std::string& linear_solver_name, const std::string& hessian_evaluation_method, bool use_trust_region,
         bool scale_residuals);

   Iterate generate_initial_point(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
   void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   void compute_progress_measures(const Problem& problem, Iterate& iterate) override;
   int get_hessian_evaluation_count() const override;

private:
   /* barrier parameter */
   double barrier_parameter;
   std::unique_ptr<HessianEvaluation> hessian_evaluation;
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that solves the subproblem */
   /* constants */
   InteriorPointParameters parameters;

   /* data structures */
   std::set<size_t> lower_bounded_variables; /* indices of the lower-bounded variables */
   std::set<size_t> upper_bounded_variables; /* indices of the upper-bounded variables */

   bool force_symbolic_factorization{true};
   double inertia_hessian;
   double inertia_hessian_last_;
   double inertia_constraints;
   double default_multiplier_;
   size_t iteration;
   size_t number_factorizations_;

   /* preallocated vectors */
   std::vector<double> rhs;
   std::vector<double> lower_delta_z;
   std::vector<double> upper_delta_z;

   void update_barrier_parameter(Iterate& current_iterate);
   void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) override;
   void factorize(COOMatrix& kkt_matrix, FunctionType problem_type);
   double compute_barrier_directional_derivative(const std::vector<double>& solution);
   double evaluate_barrier_function(const Problem& problem, Iterate& iterate);
   double compute_primal_length(const Iterate& current_iterate, const std::vector<double>& ipm_solution, double tau);
   double compute_dual_length(const Iterate& current_iterate, double tau);
   COOMatrix assemble_kkt_matrix(const Problem& problem, Iterate& current_iterate);
   void modify_inertia(COOMatrix& kkt_matrix, size_t size_first_block, size_t size_second_block, FunctionType problem_type);
   void generate_kkt_rhs(const Iterate& current_iterate);
   void compute_lower_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution);
   void compute_upper_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution);
   Direction generate_direction(const Problem& problem, const Iterate& current_iterate, std::vector<double>& solution_IPM);

   double compute_constraint_violation(const Problem& problem, const Iterate& iterate) const override;
   double compute_KKT_error_scaling(Iterate& current_iterate) const;
   double compute_predicted_reduction(const Direction& direction, double step_length) const;
   double compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds) const;
};

#endif // IPM_H
