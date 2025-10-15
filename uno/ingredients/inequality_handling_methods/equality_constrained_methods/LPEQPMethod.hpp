// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPEQPMETHOD_H
#define UNO_LPEQPMETHOD_H

#include <memory>
#include "ingredients/hessian_models/ZeroHessian.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/regularization_strategies/NoRegularization.hpp"
#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   class LPEQPMethod : public InequalityHandlingMethod {
   public:
      explicit LPEQPMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         Direction& direction, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      [[nodiscard]] double proximal_coefficient() const override;

      // matrix computations
      [[nodiscard]] EvaluationSpace& get_evaluation_space() const override;
      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

      // progress measures
      void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const OptimizationProblem& problem,
         const Iterate& iterate, const Vector<double>& primal_direction, double step_length) const override;

      void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) override;

      void set_initial_point(const Vector<double>& initial_point) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      Vector<double> initial_point{};
      const bool enforce_linear_constraints_at_initial_iterate;
      ZeroHessian LP_hessian_model{};
      NoRegularization<double> LP_regularization_strategy{};
      Direction LP_direction{};
      // pointers to allow polymorphism
      const std::unique_ptr<LPSolver> LP_solver;
      const std::unique_ptr<QPSolver> QP_solver;
      const double activity_tolerance;

      void solve_LP(Statistics& statistics, Subproblem& subproblem, const Multipliers& current_multipliers,
         const WarmstartInformation& warmstart_information);
      void solve_EQP(Statistics& statistics, Subproblem& subproblem, const Multipliers& current_multipliers,
         Direction& direction, const WarmstartInformation& warmstart_information);
   };
} // namespace

#endif // UNO_LPEQPMETHOD_H