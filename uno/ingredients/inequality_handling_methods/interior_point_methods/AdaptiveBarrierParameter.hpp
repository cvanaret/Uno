// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ADAPTIVEBARRIERPARAMETER_H
#define UNO_ADAPTIVEBARRIERPARAMETER_H

#include <algorithm>
#include <utility>
#include <vector>
#include "options/Options.hpp"

namespace uno {
   // what the interior-point method should do with the barrier parameter at the current iterate
   enum class BarrierUpdate {
      ORACLE,      // free mode: the barrier parameter is computed by the quality-function oracle
      RESET_FIXED, // switch from the free to the fixed mode: mu = monotone_init_factor * average complementarity
      MONOTONE     // fixed mode: monotone (Fiacco-McCormick) update
   };

   // mode logic of IPOPT's AdaptiveMuUpdate (Nocedal, Waechter & Waltz, SIOPT 2009) with the obj-constr-filter
   // globalization: the free mode is kept as long as the iterates make sufficient progress with respect to a filter of
   // (infeasibility, objective) pairs, which do not depend on the barrier parameter. Otherwise, the fixed mode (monotone
   // update) takes over until the iterates make sufficient progress again
   class AdaptiveBarrierParameter {
   public:
      explicit AdaptiveBarrierParameter(const Options& options):
            sigma_min(options.get_double("adaptive_mu_sigma_min")),
            sigma_max(options.get_double("adaptive_mu_sigma_max")),
            mu_min(options.get_double("adaptive_mu_min")),
            mu_max_factor(options.get_double("adaptive_mu_max_factor")),
            monotone_init_factor(options.get_double("adaptive_mu_monotone_init_factor")),
            filter_margin_factor(options.get_double("adaptive_mu_filter_margin_factor")),
            filter_max_margin(options.get_double("adaptive_mu_filter_max_margin")) {
      }

      // error: overall optimality error of the current iterate (determines the filter margin)
      [[nodiscard]] BarrierUpdate update_mode(double infeasibility, double objective, double error, double average_complementarity) {
         // the upper bound on mu is set at the first iterate
         if (this->mu_max < 0.) {
            this->mu_max = this->mu_max_factor * average_complementarity;
         }
         const bool sufficient_progress = this->makes_sufficient_progress(infeasibility, objective, error);
         if (sufficient_progress) {
            this->filter.emplace_back(infeasibility, objective);
            this->free_mode = true;
            return BarrierUpdate::ORACLE;
         }
         if (this->free_mode) {
            this->free_mode = false;
            return BarrierUpdate::RESET_FIXED;
         }
         return BarrierUpdate::MONOTONE;
      }

      [[nodiscard]] bool in_free_mode() const { return this->free_mode; }
      void leave_free_mode() { this->free_mode = false; }

      // after the restoration phase: restart in free mode with an empty filter
      void reset() {
         this->free_mode = true;
         this->filter.clear();
      }

      [[nodiscard]] double clamp(double barrier_parameter) const {
         return std::min(this->mu_max, std::max(this->mu_min, barrier_parameter));
      }

      [[nodiscard]] double fixed_mode_barrier_parameter(double average_complementarity) const {
         return this->clamp(this->monotone_init_factor * average_complementarity);
      }

      // search interval of the centering parameter sigma (mu = sigma * average complementarity)
      [[nodiscard]] std::pair<double, double> sigma_interval(double average_complementarity) const {
         const double sigma_lo = std::max(this->sigma_min, this->mu_min / average_complementarity);
         const double sigma_up = std::min(this->sigma_max, this->mu_max / average_complementarity);
         return {sigma_lo, std::max(sigma_lo, sigma_up)};
      }

   protected:
      const double sigma_min, sigma_max;
      const double mu_min, mu_max_factor;
      const double monotone_init_factor;
      const double filter_margin_factor, filter_max_margin;
      double mu_max{-1.};
      bool free_mode{true};
      std::vector<std::pair<double, double>> filter{};

      // (infeasibility, objective) is not dominated by a filter entry, up to a margin
      [[nodiscard]] bool makes_sufficient_progress(double infeasibility, double objective, double error) const {
         const double margin = this->filter_margin_factor * std::min(this->filter_max_margin, error);
         return std::none_of(this->filter.cbegin(), this->filter.cend(), [&](const auto& entry) {
            return entry.first <= infeasibility + margin && entry.second <= objective + margin;
         });
      }
   };
} // namespace

#endif // UNO_ADAPTIVEBARRIERPARAMETER_H
