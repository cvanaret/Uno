#ifndef IPM_H
#define IPM_H

#include <exception>
#include "Subproblem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "HessianModel.hpp"
#include "AugmentedSystem.hpp"

struct InteriorPointParameters {
   double tau_min;
   double k_sigma;
   double smax;
   double k_mu;
   double theta_mu;
   double k_epsilon;
};

class InteriorPoint : public Subproblem {
public:
   InteriorPoint(const Problem& problem, size_t max_number_variables, size_t number_constraints, const std::string& hessian_model,
         const std::string& linear_solver_name, const std::string& sparse_format, double initial_barrier_parameter, double default_multiplier,
         double tolerance, bool use_trust_region);
   ~InteriorPoint() override = default;

   void set_initial_point(const std::vector<double>& initial_point) override;
   void set_constraints(const Problem& problem, Iterate& iterate);
   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;
   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void add_variable(size_t i, double current_value, const Range& bounds, double objective_term, size_t j, double jacobian_term) override;
   void remove_variable(size_t i, size_t j) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   void compute_progress_measures(const Problem& problem, Iterate& iterate) override;
   void register_accepted_iterate(Iterate& iterate) override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

private:
   AugmentedSystem augmented_system;
   double barrier_parameter;
   const double tolerance;
   const std::unique_ptr<HessianModel> hessian_model;
   const std::unique_ptr<LinearSolver> linear_solver;
   const InteriorPointParameters parameters;

   // data structures
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   double default_multiplier;
   size_t iteration{0};

   // local copy of primal iterate and bound multipliers
   std::vector<double> primal_iterate;
   std::vector<double> lower_bound_multipliers;
   std::vector<double> upper_bound_multipliers;

   // preallocated vectors
   std::vector<double> barrier_constraints;
   std::vector<double> lower_delta_z;
   std::vector<double> upper_delta_z;

   void update_barrier_parameter(const Iterate& current_iterate);
   void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) override;
   double compute_barrier_directional_derivative(const std::vector<double>& solution);
   double evaluate_barrier_function(const Problem& problem, Iterate& iterate);
   double primal_fraction_to_boundary(const std::vector<double>& ipm_solution, double tau);
   double dual_fraction_to_boundary(double tau);
   void assemble_augmented_system(const Problem& problem);
   void assemble_augmented_matrix();
   void generate_augmented_rhs();
   void compute_lower_bound_dual_direction(const std::vector<double>& solution);
   void compute_upper_bound_dual_direction(const std::vector<double>& solution);
   void generate_direction(const Problem& problem, const Iterate& current_iterate);
   [[nodiscard]] double compute_KKT_error_scaling(const Iterate& current_iterate) const;
   [[nodiscard]] double compute_central_complementarity_error(const Iterate& iterate) const;
   void set_current_iterate(const Iterate& iterate);
   void print_soc_iteration(const Direction& direction_soc) const;
   void print_solution(const Problem& problem, double primal_step_length, double dual_step_length) const;
};

#endif // IPM_H
