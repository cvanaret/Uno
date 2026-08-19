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
   class Evaluations;
   class GlobalizationStrategy;
   class Iterate;
   class l1RelaxedProblem;
   class OptimizationProblem;
   class Options;
   class Statistics;
   class SolverWorkspace;
   class Subproblem;
   template <typename T>
   class Vector;
   class WarmstartInformation;
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options);
      virtual ~InequalityHandlingMethod() = default;

      virtual void generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const = 0;
      virtual void initialize_statistics(Statistics& statistics) = 0;
      [[nodiscard]] virtual bool update_parameterization(Statistics& statistics, const Iterate& current_iterate) = 0;
      [[nodiscard]] virtual const Direction& solve(Statistics& statistics, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) = 0;

      virtual void initialize_feasibility_problem(Iterate& current_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate, Evaluations& evaluations) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      [[nodiscard]] virtual bool has_second_order_corrections() const = 0;
      [[nodiscard]] virtual const Direction& compute_second_order_correction(const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) = 0;
      virtual void compute_least_squares_multipliers(Iterate& iterate, Evaluations& evaluations) = 0;

      virtual void evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const = 0;
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) const = 0;
      virtual void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      const OptimizationProblem& problem;
      const Norm progress_norm;

      void evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate, Evaluations& evaluations) const;
      bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const SolverWorkspace& solver_workspace, const Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, Evaluations& current_evaluations, Evaluations& trial_evaluations) const;
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H