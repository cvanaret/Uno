#ifndef UNO_ITERATE_H
#define UNO_ITERATE_H

#include <vector>
#include "Constraint.hpp"
#include "linear_algebra/SparseVector.hpp"

class Model;

struct Errors {
   double constraints;
   double stationarity;
   double complementarity;
};

struct ProgressMeasures {
   double infeasibility;
   double objective;
};

struct Evaluations {
   double objective{std::numeric_limits<double>::infinity()}; /*!< Objective value */
   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
   std::vector<SparseVector<double>> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

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
   std::vector<double> x; /*!< \f$\mathbb{R}^n\f$ primal variables */
   Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */

   // evaluations
   Evaluations original_evaluations;
   // lazy evaluation flags
   bool is_objective_computed{false};
   bool are_constraints_computed{false};
   bool is_objective_gradient_computed{false}; /*!< Flag that indicates if the objective gradient has already been computed */
   bool is_constraint_jacobian_computed{false}; /*!< Flag that indicates if the constraint Jacobian has already been computed */

   std::vector<double> lagrangian_gradient;

   // residuals of nonlinear functions
   Errors nonlinear_errors{0., 0., 0.};
   ProgressMeasures nonlinear_progress{0., 0.};

   static size_t number_eval_objective;
   static size_t number_eval_constraints;
   static size_t number_eval_jacobian;

   void evaluate_objective(const Model& model);
   void evaluate_constraints(const Model& model);
   void evaluate_objective_gradient(const Model& model);
   void evaluate_constraint_jacobian(const Model& model);
   void evaluate_lagrangian_gradient(const Model& model, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers);

   void set_number_variables(size_t number_variables);
   void reset_evaluations();

   friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
};

#endif // UNO_ITERATE_H
