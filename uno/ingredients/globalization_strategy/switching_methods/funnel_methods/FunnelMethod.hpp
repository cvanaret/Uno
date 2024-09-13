// Copyright (c) 2024 David Kiessling, Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FUNNELSTRATEGY_H
#define UNO_FUNNELSTRATEGY_H

#include "../SwitchingMethod.hpp"
#include "Funnel.hpp"
#include "tools/Infinity.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for funnel strategy
 *
 *  Set of constants to control the funnel strategy
 */
namespace uno {
   struct FunnelMethodParameters {
      double initial_upper_bound;
      double infeasibility_factor;
      double beta; /*!< Margin around funnel */
      double gamma; /*!< For acceptability wrt current point. Margin around objective value */
   };

   /*! \class FunnelMethod
    * \brief Step acceptance strategy based on a funnel
    *
    *  Strategy that accepts or declines a trial step
    */
   class FunnelMethod: public SwitchingMethod {
   public:
      explicit FunnelMethod(const Options& options);
      ~FunnelMethod() override = default;

      void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress,
            const ProgressMeasures& trial_progress) const override;
      void reset() override;
      void notify_switch_to_feasibility(const ProgressMeasures& current_progress_measures) override;
      void notify_switch_to_optimality(const ProgressMeasures& current_progress_measures) override;

      [[nodiscard]] bool acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
            double trial_objective) const;
      [[nodiscard]] bool infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const;
      [[nodiscard]] bool objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const;

   protected:
      Funnel funnel;
      const FunnelMethodParameters parameters; /*!< Set of constants */
      const bool require_acceptance_wrt_current_iterate;

      [[nodiscard]] bool is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) override;
      [[nodiscard]] double compute_actual_objective_reduction(double current_optimality_measure, double trial_optimality_measure);
      void set_statistics(Statistics& statistics) const override;
   };
} // namespace

#endif // UNO_FUNNELSTRATEGY_H
