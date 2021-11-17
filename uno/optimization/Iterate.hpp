#ifndef ITERATE_H
#define ITERATE_H

#include <ostream>
#include <vector>
#include "Problem.hpp"
#include "Constraint.hpp"

struct Errors {
   double constraints;
   double KKT;
   double FJ;
   double complementarity;
};

struct ProgressMeasures {
   double infeasibility;
   double objective;
};

class Iterate {
public:
   Iterate(size_t number_variables, size_t number_constraints);

   std::vector<double> x; /*!< \f$\mathbb{R}^n\f$ primal variables */
   Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */

   // functions
   double objective{std::numeric_limits<double>::infinity()}; /*!< Objective value */
   bool is_objective_computed{false};
   static size_t number_eval_objective;

   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   bool are_constraints_computed{false};
   static size_t number_eval_constraints;

   SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
   bool is_objective_gradient_computed{false}; /*!< Flag that indicates if the objective gradient has already been computed */

   std::vector<SparseVector<double>> constraints_jacobian; /*!< Sparse Jacobian of the constraints */
   bool is_constraints_jacobian_computed{false}; /*!< Flag that indicates if the constraint Jacobian has already been computed */
   static size_t number_eval_jacobian;

   std::vector<double> lagrangian_gradient;

   // residuals
   Errors errors{0., 0., 0., 0.};
   ProgressMeasures progress{0., 0.};

   void evaluate_objective(const Problem& problem);
   void evaluate_constraints(const Problem& problem);
   void evaluate_objective_gradient(const Problem& problem);
   void evaluate_constraints_jacobian(const Problem& problem);
   void evaluate_lagrangian_gradient(const Problem& problem, double objective_multiplier, const Multipliers& constraint_multipliers);

   void change_number_variables(size_t number_variables);

   friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
};

#endif // ITERATE_H
