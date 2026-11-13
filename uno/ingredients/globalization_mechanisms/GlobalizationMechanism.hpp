// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISM_H
#define UNO_GLOBALIZATIONMECHANISM_H

#include <memory>
#include <string>
#include <vector>
#include "linear_algebra/Norm.hpp"
#include "optimization/SolutionStatus.hpp"
#include "options/Options.hpp"

namespace uno {
   // forward declarations
   class ConstraintRelaxationStrategy;
   class Direction;
   class EvaluationCache;
   class Evaluations;
   class Iterate;
   class Model;
   class Statistics;
   class UserCallbacks;
   class WarmstartInformation;

   class GlobalizationMechanism {
   public:
      GlobalizationMechanism(const Model& model, bool use_trust_region, const Options& options,
         std::vector<OptionOverride>& option_overrides);
      virtual ~GlobalizationMechanism();

      virtual void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) = 0;
      virtual void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;
      virtual void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] size_t get_number_subproblems_solved() const;

      // stagnation: the accepted trial iterate differs from the current iterate by round-off only (zero or tiny
      // primal-dual step) and the parameterization (e.g. the barrier parameter) was not updated. This happens when the
      // tolerances cannot be reached (e.g. zero tolerances): the algorithm then repeats the same iterate forever. After
      // max_stagnant_iterations consecutive such iterations, the trial iterate is terminal, whatever the tolerances
      void detect_stagnation(const Model& model, const Iterate& current_iterate, Iterate& trial_iterate,
         Evaluations& trial_evaluations);
      static constexpr size_t max_stagnant_iterations{2};
      // a component v changes by round-off only if |Δv| <= stagnation_roundoff_factor * eps * (1 + |v|) (IPOPT's tiny
      // step criterion)
      static constexpr double stagnation_roundoff_factor{10.};

   protected:
      const std::unique_ptr<ConstraintRelaxationStrategy> constraint_relaxation_strategy{};
      const double primal_tolerance;
      const double tight_unscaled_primal_infeasibility_tolerance;
      const Norm residual_norm;
      size_t number_stagnant_iterations{0};
      size_t stagnation_number_parameterization_updates{0};

      // status of an iterate at which no progress can be made: acceptable status if the iterate satisfies the acceptable
      // conditions, otherwise FEASIBLE_SMALL_STEP if the iterate is feasible wrt the scaled primal tolerance or the tight
      // unscaled primal tolerance (with a zero primal tolerance, a round-off infeasibility does not make the iterate
      // infeasible), INFEASIBLE_SMALL_STEP otherwise
      [[nodiscard]] SolutionStatus small_step_status(const Model& model, const Iterate& iterate, Evaluations& evaluations) const;

      static void assemble_trial_iterate(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double primal_step_length, double constraint_dual_step_length, double bound_dual_step_length);
      static void set_primal_statistics(Statistics& statistics, const Model& model, const Iterate& iterate,
         const Evaluations& evaluations);
      static void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate);
   };
} // namespace

#endif // UNO_GLOBALIZATIONMECHANISM_H