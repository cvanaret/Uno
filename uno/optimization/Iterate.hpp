// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ITERATE_H
#define UNO_ITERATE_H

#include "SolutionStatus.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/DualResiduals.hpp"

namespace uno {
   class Iterate {
   public:
      size_t number_variables;
      size_t number_constraints;
      Vector<double> primals;
      Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */
      double objective_multiplier{1.};

      // primal-dual residuals
      double primal_feasibility{INF<double>};
      DualResiduals residuals;

      // measures of progress (infeasibility, objective, auxiliary)
      ProgressMeasures progress{INF<double>, {}, INF<double>};

      // status
      SolutionStatus status{SolutionStatus::NOT_OPTIMAL};

      // member functions
      Iterate(size_t number_variables, size_t number_constraints);
      Iterate(const Iterate& other) = default;
      Iterate(Iterate&& other) = default;
      Iterate& operator=(Iterate&& other) = default;

      void set_number_variables(size_t number_variables);

      friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
   };
} // namespace

#endif // UNO_ITERATE_H