// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ITERATE_H
#define UNO_ITERATE_H

#include <vector>
#include "ingredients/globalization_strategy/ProgressMeasures.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "optimization/Model.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/PrimalDualResiduals.hpp"
#include "tools/Infinity.hpp"

struct Evaluations {
   double objective{INF<double>}; /*!< Objective value */
   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
   RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

   Evaluations(size_t max_number_variables, size_t max_number_constraints):
         constraints(max_number_constraints),
         objective_gradient(max_number_variables),
         constraint_jacobian(max_number_constraints) {
      for (auto& constraint_gradient: this->constraint_jacobian) {
         constraint_gradient.reserve(max_number_variables);
      }
   }
};

class Iterate {
public:
   Iterate(size_t max_number_variables, size_t max_number_constraints);

   size_t number_variables;
   size_t number_constraints;
   std::vector<double> primals; /*!< \f$\mathbb{R}^n\f$ primal variables */
   Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */

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
   PrimalDualResiduals residuals{};
   LagrangianGradient<double> lagrangian_gradient;

   // measures of progress (infeasibility, optimality, auxiliary)
   ProgressMeasures progress{INF<double>, {}, INF<double>};

   // status
   TerminationStatus status{TerminationStatus::NOT_OPTIMAL};

   void evaluate_objective(const Model& model);
   void evaluate_constraints(const Model& model);
   void evaluate_objective_gradient(const Model& model);
   void evaluate_constraint_jacobian(const Model& model);

   void set_number_variables(size_t number_variables);

   friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
};

#endif // UNO_ITERATE_H
