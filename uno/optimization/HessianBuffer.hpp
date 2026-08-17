// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANBUFFER_H
#define UNO_HESSIANBUFFER_H

#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Iterate;
   class Statistics;
   class Subproblem;

   class HessianBuffer {
   public:
      HessianBuffer() = default;
      ~HessianBuffer() = default;

      void initialize(size_t number_hessian_nonzeros);
      void evaluate_hessian(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         View<double> hessian_values);

   protected:
      Vector<double> hessian_buffer{};
      bool first_evaluation{true};
   };
} // namespace

#endif // UNO_HESSIANBUFFER_H