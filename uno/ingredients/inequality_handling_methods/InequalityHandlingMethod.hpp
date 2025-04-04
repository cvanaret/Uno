// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <string>

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationSpace;
   class HessianModel;
   class Iterate;
   class l1RelaxedProblem;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   template <typename ElementType>
   class RegularizationStrategy;
   class Statistics;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod() = default;
      virtual ~InequalityHandlingMethod() = default;

      virtual void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) = 0;
      virtual void initialize_statistics(Statistics& statistics, const Options& options) = 0;
      virtual void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) = 0;
      virtual void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         Direction& direction, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) = 0;

      virtual void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
      virtual void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      // matrix computations
      [[nodiscard]] virtual EvaluationSpace& get_evaluation_space() const = 0;
      virtual void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) = 0;
      virtual void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      virtual void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      [[nodiscard]] virtual double compute_hessian_quadratic_product(const Vector<double>& vector) const = 0;

      // progress measures
      virtual void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) = 0;
      [[nodiscard]] virtual double compute_predicted_auxiliary_reduction_model(const OptimizationProblem& problem,
         const Iterate& iterate, const Vector<double>& primal_direction, double step_length) const = 0;

      virtual void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) = 0;

      virtual void set_initial_point(const Vector<double>& initial_point) = 0;
      virtual void postprocess_iterate(const OptimizationProblem& model, Vector<double>& primals, Multipliers& multipliers) = 0;

      size_t number_subproblems_solved{0};
      // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
      bool subproblem_definition_changed{false};

      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H