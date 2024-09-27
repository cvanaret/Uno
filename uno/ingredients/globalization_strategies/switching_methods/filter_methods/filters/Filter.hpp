// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTER_H
#define UNO_FILTER_H

#include <vector>
#include <memory>
#include "tools/Infinity.hpp"

namespace uno {
   // forward declaration
   class Options;

   struct FilterParameters {
      double beta; /*!< Margin around filter */
      double gamma; /*!< Margin around filter (sloping margin) */
   };

   class Filter {
   public:
      explicit Filter(const Options& options);
      virtual ~Filter() = default;

      void reset();
      [[nodiscard]] double get_smallest_infeasibility() const;
      void set_infeasibility_upper_bound(double new_upper_bound);

      virtual void add(double current_infeasibility, double current_objective);
      [[nodiscard]] virtual bool acceptable(double trial_infeasibility, double trial_objective);
      [[nodiscard]] virtual bool acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
            double trial_objective) const;
      [[nodiscard]] virtual double compute_actual_objective_reduction(double current_objective, double current_infeasibility, double trial_objective);

      [[nodiscard]] bool infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const;
      [[nodiscard]] bool objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const;

      friend std::ostream& operator<<(std::ostream& stream, Filter& filter);

   protected:
      const size_t capacity; /*!< Max filter size */
      std::vector<double> infeasibility{}; // infeasibility increases
      std::vector<double> objective{}; // objective decreases
      double infeasibility_upper_bound{INF<double>}; /*!< Upper bound on infeasibility measure */
      size_t number_entries{0};
      const FilterParameters parameters; /*!< Set of parameters */

      [[nodiscard]] bool is_empty() const;
      [[nodiscard]] bool acceptable_wrt_upper_bound(double trial_infeasibility) const;
      void left_shift(size_t start, size_t shift_size);
      void right_shift(size_t start, size_t shift_size);
   };
} // namespace

#endif // UNO_FILTER_H
