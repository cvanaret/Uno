// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "ingredients/subproblem_solvers/InequalityConstrainedSolver.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Multipliers;

   class InequalityConstrainedMethod : public InequalityHandlingMethod {
   public:
      explicit InequalityConstrainedMethod(const Options& options);
      ~InequalityConstrainedMethod() override = default;

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(Iterate& initial_iterate) override;
      void solve(Statistics& statistics, Iterate& current_iterate, Direction& direction, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      [[nodiscard]] double proximal_coefficient() const override;

      // matrix computations
      [[nodiscard]] EvaluationSpace& get_evaluation_space() const override;
      void evaluate_constraint_jacobian(Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

      // acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         UserCallbacks& user_callbacks) override;

      void postprocess_iterate(Iterate& iterate) override;

      void set_initial_point(const Vector<double>& point) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const OptimizationProblem* problem{};
      // pointer to allow polymorphism
      std::unique_ptr<InequalityConstrainedSolver> solver{};
      Vector<double> initial_point{};
      const Options& options; // copy of the options for delayed allocation of solver

      static void compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H