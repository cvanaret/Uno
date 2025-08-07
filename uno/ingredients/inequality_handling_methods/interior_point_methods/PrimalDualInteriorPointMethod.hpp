// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTMETHOD_H
#define UNO_PRIMALDUALINTERIORPOINTMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "InteriorPointParameters.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "BarrierParameterUpdateStrategy.hpp"

namespace uno {
   // forward references
   class DualResiduals;
   class Subproblem;

   class PrimalDualInteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit PrimalDualInteriorPointMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, Direction& direction,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& constraint_index) override;
      [[nodiscard]] double proximal_coefficient() const override;

      // matrix computations
      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

      // progress measures
      void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const OptimizationProblem& problem, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const override;

      void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) override;

      void set_initial_point(const Vector<double>& point) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver;
      BarrierParameterUpdateStrategy barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const double default_multiplier;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool solving_feasibility_problem{false};
      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      void update_barrier_parameter(const PrimalDualInteriorPointProblem& barrier_problem, const Iterate& current_iterate,
         const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTMETHOD_H