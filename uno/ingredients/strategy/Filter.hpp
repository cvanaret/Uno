#ifndef FILTER_H
#define FILTER_H

#include <ostream>
#include <vector>
#include <memory>
#include "tools/Options.hpp"

struct FilterConstants {
   double Beta; /*!< Margin around filter */
   double Gamma; /*!< Margin around filter (sloping margin) */
};

/*! \class Filter
 * \brief Filter
 *
 *  Filter
 */
class Filter {
public:
   explicit Filter(FilterConstants& constants);
   virtual ~Filter() = default;

   double upper_bound{std::numeric_limits<double>::infinity()}; /*!< Upper bound on constraint violation */
   const size_t max_size{50}; /*!< Max filter size */
   const FilterConstants constants; /*!< Set of constants */

   void reset();
   virtual void add(double infeasibility_measure, double optimality_measure);
   virtual bool accept(double infeasibility_measure, double optimality_measure);
   virtual bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure);
   virtual double compute_actual_reduction(double current_objective, double current_residual, double trial_objective);

   friend std::ostream& operator<<(std::ostream& stream, Filter& filter);

protected:
   std::vector<double> infeasibility;
   std::vector<double> optimality;
   size_t number_entries{0};

   void left_shift(size_t start, size_t shift_size);
   void right_shift(size_t start, size_t shift_size);
};

/*! \class NonmonotoneFilter
 * \brief Non-monotonic filter
 *
 *  Non-monotonic filter
 */
class NonmonotoneFilter : public Filter {
public:
   NonmonotoneFilter(FilterConstants& constants, size_t number_dominated_entries);

   void add(double infeasibility_measure, double optimality_measure) override;
   bool accept(double infeasibility_measure, double optimality_measure) override;
   bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure) override;
   double compute_actual_reduction(double current_objective, double current_residual, double trial_objective) override;

protected:
   const size_t max_number_dominated_entries; /*!< Memory of filter */

   size_t compute_number_dominated_entries(double infeasibility_measure, double optimality_measure);
};

class FilterFactory {
public:
   static std::unique_ptr<Filter> create(const Options& options);
};

#endif // FILTER_H
