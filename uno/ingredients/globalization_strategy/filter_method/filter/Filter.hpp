// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTER_H
#define UNO_FILTER_H

#include <vector>
#include <memory>
#include "tools/Options.hpp"
#include "tools/Infinity.hpp"

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
   [[nodiscard]] double get_infeasibility_upper_bound() const;
   void set_infeasibility_upper_bound(double new_upper_bound);

   virtual void add(double infeasibility_measure, double objective_measure);
   virtual bool acceptable(double infeasibility_measure, double objective_measure);
   virtual bool acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_objective_measure, double trial_infeasibility_measure,
         double trial_objective_measure);
   virtual double compute_actual_reduction(double current_objective_measure, double current_infeasibility_measure, double trial_objective_measure);

   friend std::ostream& operator<<(std::ostream& stream, Filter& filter);

protected:
   const size_t capacity; /*!< Max filter size */
   std::vector<double> infeasibility{}; // infeasibility increases
   std::vector<double> objective{}; // objective decreases
   double infeasibility_upper_bound{INF<double>}; /*!< Upper bound on infeasibility measure */
   size_t number_entries{0};
   const FilterParameters parameters; /*!< Set of parameters */

   [[nodiscard]] bool is_empty() const;
   [[nodiscard]] bool acceptable_wrt_upper_bound(double infeasibility_measure) const;
   void left_shift(size_t start, size_t shift_size);
   void right_shift(size_t start, size_t shift_size);
};

#endif // UNO_FILTER_H
