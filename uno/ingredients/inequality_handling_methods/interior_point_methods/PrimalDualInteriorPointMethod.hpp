// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTMETHOD_H
#define UNO_PRIMALDUALINTERIORPOINTMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "BarrierParameterUpdateStrategy.hpp"
#include "InteriorPointParameters.hpp"
#include "PrimalDualInteriorPointProblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   // forward references
   class DualResiduals;
   class Subproblem;

   class PrimalDualInteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit PrimalDualInteriorPointMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(Iterate& initial_iterate) override;
      void solve(Statistics& statistics, Iterate& current_iterate, Direction& direction, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& constraint_index) override;
      [[nodiscard]] double proximal_coefficient() const override;

      [[nodiscard]] ProgressMeasures compute_predicted_reductions(HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius, Iterate& current_iterate,
         const Direction& direction, double step_length) const override;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy,
         double trust_region_radius, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, UserCallbacks& user_callbacks) override;

      // matrix computations
      void evaluate_constraint_jacobian(Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const override;

      void postprocess_iterate(Iterate& iterate) override;

      void set_initial_point(const Vector<double>& point) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const OptimizationProblem* problem{};
      std::unique_ptr<PrimalDualInteriorPointProblem> barrier_problem{};
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      BarrierParameterUpdateStrategy barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const double default_multiplier;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      void update_barrier_parameter(const Iterate& current_iterate, const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
      [[nodiscard]] EvaluationSpace& get_evaluation_space() const override;
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTMETHOD_H