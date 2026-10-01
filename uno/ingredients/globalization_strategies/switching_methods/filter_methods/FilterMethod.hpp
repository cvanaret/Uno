// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTERMETHOD_H
#define UNO_FILTERMETHOD_H

#include "Filter.hpp"
#include "ingredients/globalization_strategies/switching_methods/SwitchingMethod.hpp"

namespace uno {
   // forward declaration
   class Options;

   struct FilterStrategyParameters {
      double upper_bound;
      double infeasibility_factor;
   };

   class FilterMethod: public SwitchingMethod {
   public:
      explicit FilterMethod(const Options& options);
      ~FilterMethod() override;

      void initialize(Statistics& statistics, const Iterate& initial_iterate) override;
      void reset() override;
      void avoid_cycling_back_to(const ProgressMeasures& current_progress) override;

   protected:
      Filter filter;
      const FilterStrategyParameters parameters; /*!< Set of constants */
   };
} // namespace

#endif // UNO_FILTERMETHOD_H