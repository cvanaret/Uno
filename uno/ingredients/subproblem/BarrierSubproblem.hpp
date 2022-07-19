// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_BARRIERSUBPROBLEM_H
#define UNO_BARRIERSUBPROBLEM_H

#include "Subproblem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "HessianModel.hpp"
#include "AugmentedSystem.hpp"
#include "tools/Options.hpp"

struct InteriorPointParameters {
   double tau_min;
   double k_sigma;
   double smax;
   double k_mu;
   double theta_mu;
   double k_epsilon;
   double barrier_update_fraction;
   double regularization_barrier_exponent;
};

class BarrierSubproblem : public Subproblem {
public:
   BarrierSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros, const Options& options);
   ~BarrierSubproblem() override = default;

   void set_initial_point(const std::optional<std::vector<double>>& optional_initial_point) override;
   void initialize(Statistics& statistics, const ReformulatedProblem& problem, Iterate& first_iterate) override;

   [[nodiscard]] double get_proximal_coefficient() const override;
   void set_elastic_variables(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] static double push_variable_to_interior(double variable_value, const Interval& variable_bounds);
   [[nodiscard]] Direction solve(Statistics& statistics, const ReformulatedProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] Direction compute_second_order_correction(const ReformulatedProblem& problem, Iterate& trial_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const ReformulatedProblem& problem, const Direction& direction) const override;
   [[nodiscard]] double compute_optimality_measure(const ReformulatedProblem& problem, Iterate& iterate) override;
   void postprocess_accepted_iterate(const ReformulatedProblem& problem, Iterate& iterate) override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

private:
   AugmentedSystem augmented_system;
   double barrier_parameter;
   double previous_barrier_parameter;
   const double tolerance;

   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */

   const std::unique_ptr<LinearSolver> linear_solver;
   const InteriorPointParameters parameters;
   const double default_multiplier;

   // preallocated vectors for bound multiplier displacements
   std::vector<double> lower_delta_z;
   std::vector<double> upper_delta_z;

   bool solving_feasibility_problem{false};

   void evaluate_functions(const ReformulatedProblem& problem, Iterate& current_iterate);
   void update_barrier_parameter(const ReformulatedProblem& problem, const Iterate& current_iterate);
   [[nodiscard]] static bool is_small_direction(const ReformulatedProblem& problem, const Iterate& current_iterate, const Direction& direction);
   [[nodiscard]] double compute_barrier_directional_derivative(const std::vector<double>& solution) const;
   [[nodiscard]] double evaluate_barrier_function(const ReformulatedProblem& problem, Iterate& iterate);
   [[nodiscard]] double primal_fraction_to_boundary(const ReformulatedProblem& problem, const Iterate& current_iterate, double tau);
   [[nodiscard]] double dual_fraction_to_boundary(const ReformulatedProblem& problem, const Iterate& current_iterate, double tau);
   void assemble_augmented_system(const ReformulatedProblem& problem, const Iterate& current_iterate);
   void assemble_augmented_matrix(const ReformulatedProblem& problem, const Iterate& current_iterate);
   void generate_augmented_rhs(const ReformulatedProblem& problem, const Iterate& current_iterate);
   void compute_lower_bound_dual_direction(const ReformulatedProblem& problem, const Iterate& current_iterate);
   void compute_upper_bound_dual_direction(const ReformulatedProblem& problem, const Iterate& current_iterate);
   void generate_direction(const ReformulatedProblem& problem, const Iterate& current_iterate);
   [[nodiscard]] double compute_KKT_error_scaling(const ReformulatedProblem& problem, const Iterate& current_iterate) const;
   [[nodiscard]] double compute_central_complementarity_error(const ReformulatedProblem& problem, const Iterate& iterate) const;
   void print_solution(const ReformulatedProblem& problem, double primal_step_length, double dual_step_length) const;
};

#endif // UNO_BARRIERSUBPROBLEM_H
