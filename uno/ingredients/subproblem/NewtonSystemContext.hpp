// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEMCONTEXT_H
#define UNO_NEWTONSYSTEMCONTEXT_H

#include "linear_algebra/View.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class Iterate;
   class Statistics;

   // values of the regularization blocks, as written into the linear system by the inertia correction strategy
   // (primal: δ_w on the regularized variables, dual: -δ_c)
   struct RegularizationValues {
      View<const double> primal;
      View<const double> dual;
   };

   // what the symbols of a Newton system may use when they are evaluated or applied at the current iterate.
   // regularization: values of the regularization blocks for the operator products. If absent, the regularization blocks
   // are applied with δ = 0
   struct NewtonSystemContext {
      Statistics& statistics;
      const Iterate& iterate;
      Evaluations& evaluations;
      const RegularizationValues* regularization{nullptr};
   };
} // namespace

#endif // UNO_NEWTONSYSTEMCONTEXT_H
