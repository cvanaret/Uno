// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ITERATE_H
#define UNO_ITERATE_H

#include "Evaluations.hpp"
#include "SolutionStatus.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/DualResiduals.hpp"

namespace uno {
   // forward declaration
   class Model;

   class Iterate {
   public:
      Iterate(size_t number_variables, size_t number_constraints);
      Iterate(const Iterate& other) = default;
      Iterate(Iterate&& other) = default;
      Iterate& operator=(Iterate&& other) = default;

      size_t number_variables;
      size_t number_constraints;
      Vector<double> primals;
      Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */
      double objective_multiplier{1.};

      // evaluations
      Evaluations evaluations;
      static size_t number_eval_objective;
      static size_t number_eval_constraints;
      static size_t number_eval_objective_gradient;
      static size_t number_eval_jacobian;
      // lazy evaluation flags
      bool is_objective_computed{false};
      bool are_constraints_computed{false};
      bool is_objective_gradient_computed{false}; /*!< Flag that indicates if the objective gradient has already been computed */
      bool is_constraint_jacobian_computed{false}; /*!< Flag that indicates if the constraint Jacobian has already been computed */

      // primal-dual residuals
      double primal_feasibility{INF<double>};
      DualResiduals residuals;

      // measures of progress (infeasibility, objective, auxiliary)
      ProgressMeasures progress{INF<double>, {}, INF<double>};

      // status
      SolutionStatus status{SolutionStatus::NOT_OPTIMAL};

      void evaluate_objective(const Model& model);
      void evaluate_constraints(const Model& model);
      void evaluate_objective_gradient(const Model& model);

      void set_number_variables(size_t number_variables);

      friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
   };
} // namespace

#endif // UNO_ITERATE_H