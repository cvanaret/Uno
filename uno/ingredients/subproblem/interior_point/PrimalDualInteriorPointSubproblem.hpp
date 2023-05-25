// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFEASIBLEINTERIORPOINTSUBPROBLEM_H
#define UNO_INFEASIBLEINTERIORPOINTSUBPROBLEM_H

#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem/HessianModel.hpp"
#include "tools/Options.hpp"
#include "BarrierParameterUpdateStrategy.hpp"

struct InteriorPointParameters {
   double tau_min;
   double k_sigma;
   double regularization_exponent;
   double small_direction_factor;
   double push_variable_to_interior_k1;
   double push_variable_to_interior_k2;
};

class PrimalDualInteriorPointSubproblem : public Subproblem {
public:
   PrimalDualInteriorPointSubproblem(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints,
         size_t max_number_jacobian_nonzeros, size_t max_number_hessian_nonzeros, const Options& options);
   ~PrimalDualInteriorPointSubproblem() override = default;

   void set_initial_point(const std::vector<double>& initial_point) override;
   void generate_initial_iterate(const NonlinearProblem& problem, Iterate& initial_iterate) override;

   void initialize_feasibility_problem() override;
   void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
   void exit_feasibility_problem(const NonlinearProblem& problem, Iterate& trial_iterate) override;

   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information) override;

   [[nodiscard]] std::function<double(double)> compute_predicted_optimality_reduction_model(const NonlinearProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const override;
   void set_auxiliary_measure(const NonlinearProblem& problem, Iterate& iterate) override;
   [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const NonlinearProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const override;

   void postprocess_iterate(const NonlinearProblem& problem, Iterate& iterate) override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

protected:
   SymmetricIndefiniteLinearSystem<double> augmented_system;
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   const std::unique_ptr<SymmetricIndefiniteLinearSolver<double>> linear_solver;

   BarrierParameterUpdateStrategy barrier_parameter_update_strategy;
   double previous_barrier_parameter;
   const double default_multiplier;
   const InteriorPointParameters parameters;
   const double least_square_multiplier_max_norm;
   const double damping_factor; // (Section 3.7 in IPOPT paper)

   // preallocated vectors for bound multiplier displacements
   std::vector<double> lower_delta_z{};
   std::vector<double> upper_delta_z{};

   bool solving_feasibility_problem{false};

   [[nodiscard]] double barrier_parameter() const;
   [[nodiscard]] double push_variable_to_interior(double variable_value, const Interval& variable_bounds) const;
   void evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information);
   void update_barrier_parameter(const NonlinearProblem& problem, const Iterate& current_iterate);
   [[nodiscard]] bool is_small_step(const NonlinearProblem& problem, const Iterate& current_iterate, const Direction& direction) const;
   [[nodiscard]] double evaluate_subproblem_objective() const;
   [[nodiscard]] double compute_barrier_term_directional_derivative(const NonlinearProblem& problem, const Iterate& current_iterate,
         const Direction& direction) const;
   [[nodiscard]] double primal_fraction_to_boundary(const NonlinearProblem& problem, const Iterate& current_iterate, double tau);
   [[nodiscard]] double dual_fraction_to_boundary(const NonlinearProblem& problem, const Iterate& current_iterate, double tau);
   void assemble_augmented_system(Statistics& statistics, const NonlinearProblem& problem, const Iterate& current_iterate);
   void generate_augmented_rhs(const NonlinearProblem& problem, const Iterate& current_iterate);
   void assemble_primal_dual_direction(const NonlinearProblem& problem, const Iterate& current_iterate);
   void compute_bound_dual_direction(const NonlinearProblem& problem, const Iterate& current_iterate);
   void compute_least_square_multipliers(const NonlinearProblem& problem, Iterate& iterate);
};

#endif // UNO_INFEASIBLEINTERIORPOINTSUBPROBLEM_H
