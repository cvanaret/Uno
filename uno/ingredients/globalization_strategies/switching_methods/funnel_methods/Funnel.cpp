// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Funnel.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"

namespace uno {
   Funnel::Funnel(const Options& options):
         margin(options.get_double("funnel_beta")),
         update_strategy(options.get_int("funnel_update_strategy")),
         kappa(options.get_double("funnel_kappa")) {
   }

   void Funnel::set_infeasibility_upper_bound(double new_upper_bound) {
      this->width = new_upper_bound;
   }

   double Funnel::current_width() const {
      return this->width;
   }

   /* check if the trial iterate is inside the current funnel */
   bool Funnel::acceptable(double trial_infeasibility) const {
      return (trial_infeasibility <= this->width);
   }

   /* check if the funnel sufficient decrease condition is satisfied */
   bool Funnel::sufficient_decrease_condition(double trial_infeasibility) const {
      return (trial_infeasibility <= this->margin * this->width);
   }

   void Funnel::update(double current_infeasibility, double trial_infeasibility) {
      if (this->update_strategy == 1) {
         if (trial_infeasibility <= current_infeasibility) {
            this->width = std::max(this->margin * this->width,
                  Funnel::convex_combination(current_infeasibility, trial_infeasibility, this->kappa));
         }
         else {
            DEBUG << "Trial infeasibility higher than current infeasibility" << '\n';
            this->width = this->margin * this->width;
         }
      }
      else if (this->update_strategy == 2) {
         this->width = Funnel::convex_combination(this->width, trial_infeasibility, this->kappa);
      }
      else if (this->update_strategy == 3) {
         this->width = this->margin * this->width;
      }
      else {
         throw std::runtime_error("Funnel update strategy " + std::to_string(this->update_strategy) + " is unknown");
      }
      DEBUG << "\t\tNew funnel parameter is: " << this->width << '\n';
   }

   void Funnel::update_restoration(double current_infeasibility) {
      this->width = Funnel::convex_combination(this->width, current_infeasibility, this->kappa);
      DEBUG << "\t\tNew funnel parameter is: " << this->width << '\n';
   }

   void Funnel::print() const {
      DEBUG << "Current funnel width: " << this->width << '\n';
   }

   double Funnel::convex_combination(double a, double b, double coefficient) {
      return coefficient*a + (1. - coefficient)*b;
   }
} // namespace