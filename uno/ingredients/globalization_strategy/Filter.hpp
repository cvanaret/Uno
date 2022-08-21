// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTER_H
#define UNO_FILTER_H

#include <vector>
#include <memory>
#include "tools/Options.hpp"
#include "tools/Infinity.hpp"

struct FilterConstants {
   double beta; /*!< Margin around filter */
   double gamma; /*!< Margin around filter (sloping margin) */
};

class Filter {
public:
   double upper_bound{INF<double>}; /*!< Upper bound on constraint violation */

   explicit Filter(const Options& options);
   virtual ~Filter() = default;

   void reset();
   virtual void add(double infeasibility_measure, double optimality_measure);
   virtual bool accept(double infeasibility_measure, double optimality_measure);
   virtual bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure);
   virtual double compute_actual_reduction(double current_optimality_measure, double current_infeasibility_measure, double trial_optimality_measure);

   friend std::ostream& operator<<(std::ostream& stream, Filter& filter);

protected:
   const size_t capacity; /*!< Max filter size */
   std::vector<double> infeasibility{}; // infeasibility increases
   std::vector<double> optimality{}; // optimality decreases
   size_t number_entries{0};
   const FilterConstants constants; /*!< Set of constants */

   void left_shift(size_t start, size_t shift_size);
   void right_shift(size_t start, size_t shift_size);
};

class NonmonotoneFilter : public Filter {
public:
   explicit NonmonotoneFilter(const Options& options);

   void add(double infeasibility_measure, double optimality_measure) override;
   bool accept(double infeasibility_measure, double optimality_measure) override;
   bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure) override;
   double compute_actual_reduction(double current_optimality_measure, double current_infeasibility_measure, double trial_optimality_measure) override;

protected:
   const size_t max_number_dominated_entries; /*!< Memory of filter */

   size_t compute_number_dominated_entries(double infeasibility_measure, double optimality_measure);
};

class FilterFactory {
public:
   static std::unique_ptr<Filter> create(const Options& options);
};

#endif // UNO_FILTER_H
