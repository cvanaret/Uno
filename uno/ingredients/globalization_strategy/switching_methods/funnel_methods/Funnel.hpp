// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FUNNEL_H
#define UNO_FUNNEL_H

#include "tools/Infinity.hpp"

namespace uno {
   // forward reference
   class Options;

   class Funnel {
   public:
      explicit Funnel(const Options& options);
      ~Funnel() = default;

      void set_infeasibility_upper_bound(double new_upper_bound);
      [[nodiscard]] double current_width() const;
      [[nodiscard]] bool acceptable(double trial_infeasibility) const;
      [[nodiscard]] bool sufficient_decrease_condition(double trial_infeasibility) const;
      void update(double current_infeasibility, double trial_infeasibility);
      void update_restoration(double current_infeasibility);

      void print() const;

   protected:
      double width{INF<double>};
      const double margin;
      const int update_strategy;
      const double kappa;

      [[nodiscard]] static double convex_combination(double a, double b, double coefficient);
   };
} // namespace

#endif // UNO_FUNNEL_H