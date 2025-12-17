// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"

namespace uno {
   class BacktrackingLineSearch : public GlobalizationMechanism {
   public:
      BacktrackingLineSearch(const Model& model, const Options& options);
      ~BacktrackingLineSearch() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction) override;
      void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      [[nodiscard]] std::string get_name() const override;

   private:
      const double backtracking_ratio;
      const double minimum_step_length;
      const bool scale_duals_with_step_length;

      void backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const;
      [[nodiscard]] bool terminate_with_small_step_length(Statistics& statistics, const Model& model, Iterate& trial_iterate) const;
      [[nodiscard]] double decrease_step_length(double step_length) const;
      static void check_unboundedness(const Direction& direction);

      static void set_LS_statistics(Statistics& statistics, size_t number_iterations);
   };
} // namespace

#endif // UNO_BACKTRACKINGLINESEARCH_H