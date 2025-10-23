// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <functional>
#include <string>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "linear_algebra/Norm.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationSpace;
   class GlobalizationStrategy;
   class HessianModel;
   class Iterate;
   class l1RelaxedProblem;
   class Model;
   class OptimizationProblem;
   class Options;
   template <typename ElementType>
   class InertiaCorrectionStrategy;
   class Statistics;
   class Subproblem;
   class UserCallbacks;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;
   
   class InequalityHandlingMethod {
   public:
      explicit InequalityHandlingMethod(const Options& options);
      virtual ~InequalityHandlingMethod() = default;

      virtual void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) = 0;
      virtual void initialize_statistics(Statistics& statistics, const Options& options) = 0;
      virtual void generate_initial_iterate(Iterate& initial_iterate) = 0;
      virtual void solve(Statistics& statistics, Iterate& current_iterate, Direction& direction, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) = 0;

      virtual void initialize_feasibility_problem(Iterate& current_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      // matrix computations
      [[nodiscard]] virtual EvaluationSpace& get_evaluation_space() const = 0;
      virtual void evaluate_constraint_jacobian(Iterate& iterate) = 0;
      virtual void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      virtual void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      [[nodiscard]] virtual double compute_hessian_quadratic_product(const Vector<double>& vector) const = 0;

      // progress measures
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         UserCallbacks& user_callbacks) = 0;

      virtual void postprocess_iterate(Iterate& iterate) = 0;

      virtual void set_initial_point(const Vector<double>& initial_point) = 0;

      size_t number_subproblems_solved{0};

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      const Norm progress_norm;
      // first_order_predicted_reduction is true when the predicted reduction can be taken as first-order (e.g. in line-search methods)
      const bool first_order_predicted_reduction;
      // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
      bool subproblem_definition_changed{false};

      void evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate) const;
      [[nodiscard]] double compute_predicted_infeasibility_reduction(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const;
      [[nodiscard]] std::function<double(double)> compute_predicted_objective_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const;
      [[nodiscard]] ProgressMeasures compute_predicted_reductions(const Subproblem& subproblem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, UserCallbacks& user_callbacks);
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H