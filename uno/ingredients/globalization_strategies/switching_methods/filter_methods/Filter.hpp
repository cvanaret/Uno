// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTER_H
#define UNO_FILTER_H

#include <string>
#include <vector>
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
      ~Filter() = default;

      // empty the filter and restore the upper bound set by set_infeasibility_upper_bound
      void reset();
      [[nodiscard]] double get_smallest_infeasibility() const;
      // set the upper bound and the value restored upon reset
      void set_infeasibility_upper_bound(double new_upper_bound);

      [[nodiscard]] bool acceptable_wrt_infeasibility_upper_bound(double trial_infeasibility) const;
      [[nodiscard]] bool filter_acceptable(double trial_infeasibility, double trial_objective) const;
      [[nodiscard]] bool acceptable_wrt_current_iterate(double current_infeasibility, double current_objective,
         double trial_infeasibility, double trial_objective) const;
      [[nodiscard]] bool infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const;

      void add(double current_infeasibility, double current_objective);

      friend std::ostream& operator<<(std::ostream& stream, const Filter& filter);

   protected:
      // maximum number of entries. 0: unlimited capacity (the filter grows)
      const size_t capacity;
      std::vector<double> infeasibility{}; // infeasibility increases
      std::vector<double> objective{}; // objective decreases
      double infeasibility_upper_bound{Inf}; /*!< Upper bound on infeasibility measure */
      // the upper bound is lowered when an entry is dropped from a full filter; this value is restored upon reset
      double initial_infeasibility_upper_bound{Inf};
      const FilterParameters parameters; /*!< Set of parameters */
      const bool margin_uses_entry_infeasibility;
      static constexpr size_t fixed_length_column1 = 14;
      static constexpr size_t fixed_length_column2 = 10;

      [[nodiscard]] size_t number_entries() const;
      [[nodiscard]] bool is_empty() const;
      [[nodiscard]] bool objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const;
      static void print_line(std::string& table, const std::string& infeasibility, const std::string& objective);
      [[nodiscard]] bool is_sorted() const;
   };
} // namespace

#endif // UNO_FILTER_H