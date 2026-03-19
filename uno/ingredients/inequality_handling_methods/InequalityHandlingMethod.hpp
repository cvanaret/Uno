// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <memory>
#include <string>
#include "linear_algebra/Norm.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationCache;
   class Evaluations;
   class SolverWorkspace;
   class GlobalizationStrategy;
   class HessianModel;
   class InertiaCorrectionStrategy;
   class Iterate;
   class l1RelaxedProblem;
   class OptimizationProblem;
   class Options;
   class Parameterization;
   class Statistics;
   class Subproblem;
   class UserCallbacks;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod() = default;
      virtual ~InequalityHandlingMethod() = default;

      virtual void check_problem(const OptimizationProblem& problem, bool uses_trust_region) = 0;
      virtual void initialize_statistics(Statistics& statistics) = 0;
      [[nodiscard]] virtual std::unique_ptr<OptimizationProblem> reformulate(const OptimizationProblem& problem,
         Parameterization& parameterization) = 0;
      virtual void update_parameterization(Statistics& statistics, const OptimizationProblem& problem,
         const Iterate& current_iterate, Parameterization& parameterization) = 0;

      virtual void initialize_feasibility_problem(Iterate& current_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate, Evaluations& evaluations) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const SolverWorkspace& solver_workspace, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, EvaluationCache& evaluation_cache, UserCallbacks& user_callbacks);

      size_t number_subproblems_solved{0};

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
      bool subproblem_definition_changed{false};
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H